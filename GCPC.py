import os
import csv
import argparse

def GCcount(seqUp):
    GC = 0
    lenSeq = len(seqUp)
    for i in range(lenSeq):
        if seqUp[i] == 'G' or seqUp[i] == 'C':
            GC += 1
    return GC/lenSeq

def countPAM(seqUp):
    PAMcountDic = dict()
    for i in range(len(seqUp)-3):
        if ((seqUp[i] == 'G' or seqUp[i] == 'T') and (seqUp[i+1] == 'C' or seqUp[i+1] == 'T') and seqUp[i+2] == 'T' and (seqUp[i+3] == 'A' or seqUp[i+3] == 'C' or seqUp[i+3] == 'G')) \
        or ((seqUp[i+3] == 'C' or seqUp[i+3] == 'A') and (seqUp[i+2] == 'G' or seqUp[i+2] == 'A') and seqUp[i+1] == 'A' and (seqUp[i] == 'T' or seqUp[i] == 'C' or seqUp[i] == 'G')) :
            PAMcountDic['KYTV'] = PAMcountDic.get('KYTV', 0) + 1
        if seqUp[i:i+2] == 'GG' or seqUp[i:i+2] == 'CC':
            PAMcountDic['NGG'] = PAMcountDic.get('NGG', 0) + 1
        if ((seqUp[i] == 'C' or seqUp[i] == 'T') and seqUp[i+1] == 'T') or (seqUp[i+1] == 'A' or (seqUp[i] == 'G' or seqUp[i] == 'A')):
            PAMcountDic['YTN'] = PAMcountDic.get('YTN', 0) + 1
        if ((seqUp[i] == 'G' or seqUp[i] == 'T') and (seqUp[i+1: i+3] == 'TT') and (seqUp[i+3] == 'A' or seqUp[i+3] == 'C' or seqUp[i+3] == 'G')) \
               or ((seqUp[i+3] == 'C' or seqUp[i+3] == 'A') and (seqUp[i+1: i+3] == 'AA') and (seqUp[i] == 'A' or seqUp[i] == 'C' or seqUp[i] == 'G')):
            PAMcountDic['KTTV'] = PAMcountDic.get('KTTV', 0) + 1
        if seqUp[i:i+2] == 'TT' or seqUp[i:i+2] == 'AA':
            PAMcountDic['TTN'] = PAMcountDic.get('TTN', 0) + 1
    return PAMcountDic


def mainStep(fastaDir, savePath):
    seq = ''
    f = open(savePath, 'w', newline='')
    rowNameLst = ['Strain', 'GC%', 'TTN',  'YTN', 'KYTV', 'KTTV', 'NGG']
    fCsv = csv.writer(f)
    fCsv.writerow(rowNameLst)
    print('Start Calculating...')
    for fastaName in os.listdir(fastaDir):
        if fastaName.endswith('.fasta'):
            seq, prompt = '', ''
            strain = fastaName[:-6]
            with open(os.path.join(fastaDir, fastaName)) as fhand:
                for line in fhand:
                    if not line.startswith(">"):
                        seq += line.strip()
                seqUp = seq.upper()
                GC = GCcount(seqUp)
                PAMcountDic = countPAM(seqUp)
                rowDic = PAMcountDic
                rowDic['Strain'] = strain
                rowDic['GC%'] = '{:.4f}'.format(GC)
                rowLSt = list()
                for name in rowNameLst:
                    prompt += '%s: %s\t'%(name, rowDic[name])
                    rowLSt.append(rowDic[name])
                print(prompt)
                fCsv.writerow(rowLSt)
    f.close()
    print('Complete Calculation')


parser = argparse.ArgumentParser(description="GCPC (GC content and PAM Calculation), is a tool for calculating the GC content and PAM count in different strains.\n\t\n\tThe tool support five types of PAMs, including 'TTN', 'YTN', 'KYTV', 'KTTV', and 'NGG'.\n\tThe output format is '.csv'. \n\t\n\tIf you demand other PAMs calculation, please mailto: wangzhp@shanghaitech.edu.cn", formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('-d', '--genomeDir', nargs=1, metavar='', help="The directory that store the genome files (only support '.fasta' format)", type=str)
parser.add_argument('-o', '--output', nargs=1, default='GCPC_output', metavar='', help="The output fileName (default name is 'GCPC_output)'. The output is under the genome directory.", type=str)


if __name__ == '__main__':
    args = parser.parse_args()
    genomeDir = args.genomeDir[0]
    savePath = os.path.join(genomeDir, args.output + '.csv')
    mainStep(genomeDir, savePath)
