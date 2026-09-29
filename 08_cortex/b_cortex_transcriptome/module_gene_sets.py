# Pre-specified pathway gene sets: complement cascade in seven layers (43 genes) plus the other pathways of the 12 panel proteins.
# Inputs (relative to ASD_ROOT): none
# Outputs (relative to ASD_ROOT): panel12/data/gene_modules.tsv  (module, gene, in_panel12)
# Usage: python module_gene_sets.py
# Read by gandal2022_primary.py, replication_modules.py and ../c_single_nucleus/wamsley_pseudobulk.py.
import csv
import os

ROOT = os.environ.get("ASD_ROOT", ".")
D = os.path.join(ROOT, 'panel12', 'data')
TAB, NL = chr(9), chr(10)

MODULES = {
    # complement, by cascade layer (CMP_ prefix = complement)
    'CMP_1_recognition_classical': ['C1QA', 'C1QB', 'C1QC', 'C1R', 'C1S', 'C1RL'],
    'CMP_2_recognition_lectin': ['MBL2', 'MASP1', 'MASP2', 'FCN1', 'FCN2', 'FCN3',
                                 'COLEC11'],
    'CMP_3_alternative': ['CFB', 'CFD', 'CFP'],
    'CMP_4_convertase': ['C2', 'C3', 'C4A', 'C4B', 'C5'],
    'CMP_5_MAC': ['C6', 'C7', 'C8A', 'C8B', 'C8G', 'C9'],
    'CMP_6_regulators': ['CFH', 'CFI', 'CD46', 'CD55', 'CD59', 'SERPING1', 'CLU',
                         'VTN', 'CFHR1'],
    'CMP_7_receptors': ['C3AR1', 'C5AR1', 'CR1', 'ITGAM', 'ITGB2', 'CD93', 'VSIG4'],
    # other pathways of the panel proteins
    'ECM_integrin': ['POSTN', 'COL1A1', 'COL1A2', 'COL4A1', 'COL6A1', 'COL6A2',
                     'COL6A3', 'FN1', 'LUM', 'COMP', 'TNC', 'SPARC', 'SPARCL1',
                     'MMP2', 'MMP9', 'TIMP1', 'TIMP2', 'ITGAV', 'ITGA1', 'ITGA2',
                     'ITGA5', 'ITGB1', 'ITGB5', 'DDR1', 'DDR2', 'SDC1', 'SDC4',
                     'LRP1', 'CD44'],
    'IGF_axis': ['IGF1', 'IGF2', 'IGF1R', 'IGF2R', 'IGFBP2', 'IGFBP3', 'IGFBP4',
                 'IGFBP5', 'IGFBP6', 'IGFBP7', 'IGFALS', 'INSR'],
    'Peptidase_redox': ['ANPEP', 'DPP4', 'MME', 'BCHE', 'ACHE', 'BTD', 'QSOX1',
                        'QSOX2', 'TXN', 'TXNRD1', 'PRDX1', 'PRDX6', 'GPX1',
                        'SOD1', 'SOD2'],
    'Prostaglandin': ['PTGDS', 'HPGDS', 'PTGDR', 'PTGDR2', 'PTGER2', 'PTGER3',
                      'PTGES3', 'PTGS1', 'PTGS2'],
    # cell-state context
    'Microglia_state': ['AIF1', 'CX3CR1', 'P2RY12', 'TMEM119', 'TYROBP', 'TREM2',
                        'CSF1R', 'SPI1'],
    'Synapse': ['SYP', 'SNAP25', 'DLG4', 'SYN1', 'GRIN1', 'GRIA1', 'NRXN1',
                'NLGN3', 'SHANK3', 'SYT1'],
    'BBB_transport': ['SLC2A1', 'CLDN5', 'PECAM1', 'ABCB1', 'MFSD2A', 'TFRC',
                      'LRP1', 'SLC7A5'],
}
PANEL12 = ['AHSG', 'ANPEP', 'BCHE', 'BTD', 'C1RL', 'C3', 'CLEC3B', 'IGFBP5',
           'MBL2', 'POSTN', 'PTGDS', 'QSOX1']

if __name__ == '__main__':
    os.makedirs(D, exist_ok=True)
    with open(os.path.join(D, 'gene_modules.tsv'), 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f, delimiter=TAB, lineterminator=NL)
        w.writerow(['module', 'gene', 'in_panel12'])
        for m, gs in MODULES.items():
            for g in gs:
                w.writerow([m, g, g in PANEL12])
    n_cmp = sum(len(v) for k, v in MODULES.items() if k.startswith('CMP_'))
    print('modules %d | complement genes %d' % (len(MODULES), n_cmp))
