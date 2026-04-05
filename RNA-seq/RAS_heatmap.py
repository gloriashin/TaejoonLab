import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from matplotlib.patches import Patch
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist

exclude_PC1 = ["st14-NNEa-B108T92", "st14-NNEp-B108T95"]
exclude_PC1_2 = ["st12-NBa-B108T89", "st12-NPp-B100T16",
                 "st14-NP-B100T18", "st14-NNEa-B108T92", "st14-NNEp-B108T95", "st17-NC-B100T3"]


def geneid_index(csv_file='merged_tpm_data.csv', stage='12'):
    """
    병합된 CSV 파일을 불러와서 열기
    """
    # CSV 파일 읽기 (탭으로 구분된 파일)
    df = pd.read_csv(csv_file)
    print(f"\n데이터 로드 완료: {df.shape[0]} genes x {df.shape[1]-1} samples")
    
    # Gene ID에서 'GeneID:' 제거
    df['Gene ID'] = df['Gene ID'].astype(str).str.replace('GeneID:', '', regex=False)
    
    # 결측치 처리 (0으로 채움)
    df = df.fillna(0)
    
    return df


def exclude(df, exclude_list):
    """
    exclude 리스트에 해당하는 column 없애기 
    """
    if len(exclude_list) > 0:
        # 실제로 존재하는 컬럼만 제외
        existing_cols = [col for col in exclude_list if col in df.columns]
        if existing_cols:
            mod_df = df.drop(labels=existing_cols, axis=1, errors='ignore')
            print(f"제외된 샘플: {existing_cols}")
        else:
            mod_df = df.copy()
            print("제외할 샘플이 데이터에 없습니다.")
    else:
        mod_df = df.copy()
    
    return mod_df


def analyze_ras_expression(df):
    """
    RAS 유전자 발현 분석 함수
    
    Parameters:
    -----------
    df : pandas.DataFrame
        'Gene ID' 컬럼과 샘플 컬럼들을 포함하는 데이터프레임
    
    Returns:
    --------
    processed_data : pandas.DataFrame
        필터링 및 처리된 RAS 유전자 데이터 (wide format)
    """
    
    # RAS 유전자 Entrez ID 목록
    ras_entrez_ids = ["399406", "100505423", "108712520", "399446"]
    
    # RAS 유전자 필터링
    ras_data = df[df['Gene ID'].isin(ras_entrez_ids)].copy()
    
    if ras_data.empty:
        print("경고: RAS 유전자를 찾을 수 없습니다!")
        return pd.DataFrame()
    
    # Gene ID를 Gene_Name으로 변환
    gene_name_mapping = {
        "399406": "HRAS",
        "100505423": "KRAS.L",
        "108712520": "KRAS.S",
        "399446": "NRAS"
    }
    
    ras_data['Gene_Name'] = ras_data['Gene ID'].map(gene_name_mapping)
    
    # Gene_Name을 인덱스로 설정하고 Gene ID 컬럼 제거
    ras_data = ras_data.set_index('Gene_Name')
    ras_data = ras_data.drop('Gene ID', axis=1)
    
    print(f"발견된 RAS 유전자: {list(ras_data.index)}")
    
    return ras_data


def create_heatmap(processed_data, title="RAS Paralog Expression Heatmap", 
                   figsize=(14, 7), cmap="RdBu_r"):
    """
    RAS 유전자 발현 히트맵 생성 (Log TPM 값 사용)
    
    Parameters:
    -----------
    processed_data : pandas.DataFrame
        analyze_ras_expression 함수의 출력
    title : str
        히트맵 제목
    figsize : tuple
        그림 크기
    """
    
    if processed_data.empty:
        print("데이터가 비어있어 히트맵을 생성할 수 없습니다.")
        return None
    
    # Log transformation (log2(TPM + 1))
    heatmap_data_log = np.log2(processed_data + 1)
    
    # 클러스터링을 위한 거리 계산
    # Row clustering
    if len(heatmap_data_log) > 1:
        row_linkage = hierarchy.linkage(
            pdist(heatmap_data_log, metric='euclidean'), 
            method='complete'
        )
    else:
        row_linkage = None
    
    # Column clustering
    if len(heatmap_data_log.columns) > 1:
        col_linkage = hierarchy.linkage(
            pdist(heatmap_data_log.T, metric='euclidean'), 
            method='complete'
        )
    else:
        col_linkage = None
    
    # ClusterMap 생성
    g = sns.clustermap(
        heatmap_data_log,
        cmap=cmap,
        center=heatmap_data_log.median().median(),
        row_linkage=row_linkage,
        col_linkage=col_linkage,
        cbar_kws={'label': 'log2(TPM + 1)'},
        xticklabels=True,
        yticklabels=True,
        figsize=figsize,
        dendrogram_ratio=0.15,
        cbar_pos=(0.02, 0.8, 0.03, 0.15)
    )
    
    # 제목 추가
    g.fig.suptitle(title + " (Log TPM)", fontsize=14, y=1.02)
    
    # X축 레이블 회전
    plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=45, ha='right')
    
    # 저장
    output_filename = title.replace(' ', '_') + '.png'
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"Heatmap saved as '{output_filename}'")
    
    return g


def create_heatmap_z(processed_data, title="RAS Paralog Expression Heatmap zscore", 
                   figsize=(14, 7), cmap="RdBu_r"):
    """
    RAS 유전자 발현 히트맵 생성 (Z-score 사용)
    
    Parameters:
    -----------
    processed_data : pandas.DataFrame
        analyze_ras_expression 함수의 출력
    title : str
        히트맵 제목
    figsize : tuple
        그림 크기
    """
    
    if processed_data.empty:
        print("데이터가 비어있어 히트맵을 생성할 수 없습니다.")
        return None
    
    # Log transformation (log2(TPM + 1))
    heatmap_data_log = np.log2(processed_data + 1)
    
    # Z-score 계산 (각 행에 대해)
    heatmap_data_zscore = heatmap_data_log.T  # 전치
    heatmap_data_zscore = (heatmap_data_zscore - heatmap_data_zscore.mean()) / heatmap_data_zscore.std()
    heatmap_data_zscore = heatmap_data_zscore.T  # 다시 전치
    
    # 클러스터링을 위한 거리 계산
    # Row clustering
    if len(heatmap_data_zscore) > 1:
        row_linkage = hierarchy.linkage(
            pdist(heatmap_data_zscore, metric='euclidean'), 
            method='complete'
        )
    else:
        row_linkage = None
    
    # Column clustering
    if len(heatmap_data_zscore.columns) > 1:
        col_linkage = hierarchy.linkage(
            pdist(heatmap_data_zscore.T, metric='euclidean'), 
            method='complete'
        )
    else:
        col_linkage = None
    
    # ClusterMap 생성
    g = sns.clustermap(
        heatmap_data_zscore,
        cmap=cmap,
        center=0,
        row_linkage=row_linkage,
        col_linkage=col_linkage,
        cbar_kws={'label': 'Z-score'},
        xticklabels=True,
        yticklabels=True,
        figsize=figsize,
        dendrogram_ratio=0.15,
        cbar_pos=(0.02, 0.8, 0.03, 0.15),
        vmin=-2,
        vmax=2
    )
    
    # 제목 추가 (y 위치를 1.02로 조정하여 dendogram과 겹치지 않게)
    g.fig.suptitle(title + " (Z-score)", fontsize=14, y=1.02)
    
    # X축 레이블 회전
    plt.setp(g.ax_heatmap.xaxis.get_majorticklabels(), rotation=45, ha='right')
    
    # 저장
    output_filename = title.replace(' ', '_') + '_zscore.png'
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    print(f"Heatmap saved as '{output_filename}'")
    
    return g

if __name__ == '__main__':
    stagelist = ['12', '14', '17']
    
    for num in stagelist:
        outputcsv = 'merged_tpm_data_' + num + '.csv'
        
        # 파일 존재 확인
        if not os.path.exists(outputcsv):
            print(f"경고: {outputcsv} 파일이 존재하지 않습니다. 건너뜁니다.")
            continue
        
        print(f"\n{'='*60}")
        print(f"Stage {num} 분석 시작")
        print(f"{'='*60}")
        
        # CSV 파일을 DataFrame으로 로드
        base_df = geneid_index(outputcsv, num)
        
        # Part 2 실행
        exclude_list = [[], exclude_PC1, exclude_PC1_2]
        exclude_name = ["_original", "_excludePC1", "_excludePC1_2"]
        
        for i in range(len(exclude_list)):
            print(f"\n--- {exclude_name[i]} 분석 ---")
            
            # DataFrame에서 제외할 샘플 제거
            df = exclude(base_df, exclude_list[i])
            
            # RAS 유전자 추출
            ras_df = analyze_ras_expression(df)
            
            # 히트맵 생성
            if not ras_df.empty:
                title = "RAS_Paralog_Expression_Heatmap_" + num + exclude_name[i]
                create_heatmap(ras_df, title,cmap='Reds')
                plt.close()  # 메모리 관리를 위해 figure 닫기
                
                create_heatmap_z(ras_df, title)
                plt.close()
                
            
    print("\n모든 분석이 완료되었습니다!")

