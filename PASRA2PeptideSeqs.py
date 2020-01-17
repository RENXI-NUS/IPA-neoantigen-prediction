import sys
import numpy as np
import subprocess
from Bio.Seq import Seq
import bisect

def createUniqueIntronList(csvfile, outpath):
    # Read in chromosome locations and TPM values (columns 1 and 4) from KMA output file
    chromlocs = np.loadtxt(csvfile, dtype=str, delimiter=',', skiprows=1, usecols=[1,3])
    # Only extract unique chromosomal locations
    _, indices = np.unique(chromlocs[:,0], return_index=True)
    uniquelocs = chromlocs[indices,:]
    uniquelocs = np.core.defchararray.strip(uniquelocs, '"')  # Strip "s from locations
    # Write contents of uniquelocs array to file (for reference later)
    filepath = outpath+'/'+'uniqueIntronList.txt'
    np.savetxt(filepath, uniquelocs, fmt='%s', delimiter='\t')
    # Return unique list of intron locations and corresponding TPM values
    return list(uniquelocs[:,0]),list(uniquelocs[:,1])
        
def manualTranslate(fastasequence):
    # Initialize codon table and list of peptides
    codontable = {'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*', 'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}

        # Translate full length fasta sequence all the way through
    fulllengthprotein = ''
    for i in xrange(0, len(fastasequence), 3):
        codon = fastasequence[i:i+3]
        # Account for bizarre edge cases that should really never happen
        if len(codon) != 3 or codon not in codontable:
            break
        if codontable[codon] == '*':  # Stop translating when we hit a stop codon
            break
        AA = codontable[codon]
        fulllengthprotein += AA
    return fulllengthprotein

def getSeqs(nAAs, outpath):
    outfile = open(outpath + '/peptideSeqsFASTA_' + sys.argv[3] + '.fa','a')
    headermapfile = open(outpath + '/headermap_' + sys.argv[3],'a')
    uniqueIntrons = open(sys.argv[3], 'r')
    for loc in uniqueIntrons:
        uniqueIntrons = open(sys.argv[3], 'r')
        chrom = loc.split("\t")[0]
        peptideends = loc.split("\t")[1]
        peptideends = int(peptideends)
        #print peptideends
        intron_index = loc.split("\t")[2].split('_')[3]
        intron_index = int(intron_index)
        if intron_index == 0:
            continue
        exonstarts = loc.split("\t")[3].split(',')[intron_index]
        exonstarts = int(exonstarts) + 1
        exonends = loc.split("\t")[4].split(',')[intron_index]
        exonends = int(exonends)
        exonframes = loc.split("\t")[5].split(',')[intron_index]
        exonframes = int(exonframes)
        strand = loc.split("\t")[6]
        #       intronstart = loc.split("\t")[6]
        #intronstart = int(intronstart)
        tpmval = loc.split("\t")[8] 
        #tpmval = int(tpmval)
        gene_name = loc.split("\t")[9]
        #file_name = loc.split("\t")[10]
        if strand == ('-'):
            intron_index = intron_index + 1
            element_number = len(loc.split("\t")[5].split(','))
            if element_number <= intron_index:
                continue
            exonstarts = loc.split("\t")[3].split(',')[intron_index]
            exonstarts = int(exonstarts) + 1
        #    print exonstarts
            exonends = loc.split("\t")[4].split(',')[intron_index]
            exonends = int(exonends)
        #    print exonends
            exonframes = loc.split("\t")[5].split(',')[intron_index]
            exonframes = int(exonframes)
        # Check to make sure exon starts, ends, and frames list are the same length, and if there is an error, skip this intron
        #if not len(exonstarts) == len(peptideends) == len(exonframes):
                #        continue
                # Get ORF orientation at the start of the intron
                # Handle + and - strand cases separately (need to look for different positions)
        ORForientation = 0
        #if strand == ('+'):
        if strand == ('+'):
            if exonframes == 0:
                myframe = (exonends-exonstarts) % 3
                if myframe == 0:
                    ORForientation = 0
                elif myframe == 1:
                    ORForientation = 2
                elif myframe == 2:
                    ORForientation = 1
                else: # If frame == -1 (meaning no translation takes place according to table browser)
                    continue
            elif exonframes == 1:
                myframe = (exonends-(exonstarts+2)) % 3
                if myframe == 0:
                    ORForientation = 0
                elif myframe == 1:
                    ORForientation = 2
                elif myframe == 2:
                    ORForientation = 1    
                else:
                    continue
            elif exonframes == 2:
                myframe = (exonends-(exonstarts+1)) % 3
                if myframe == 0:
                    ORForientation = 0
                elif myframe == 1:
                    ORForientation = 2
                elif myframe == 2:
                    ORForientation = 1
                else:
                    continue
            else:
                continue
        else:
                #        myframe = (exonstarts-exonframes+nAAs*3-1)) % 3
            if exonframes == 0:
                myframe = (exonends-exonstarts) % 3
                if myframe == 0:
                    ORForientation = 0
                elif myframe == 1:
                    ORForientation = 2
                elif myframe == 2:
                    ORForientation = 1
                else:
                    continue
            elif exonframes == 1:
                myframe = (exonends-2-exonstarts) % 3
                if myframe == 0:
                    ORForientation = 0
                elif myframe == 1:
                    ORForientation = 2
                elif myframe == 2:
                    ORForientation = 1
                else:
                    continue
            elif exonframes == 2:
                myframe = (exonends-1-exonstarts) % 3
                if myframe == 0:
                    ORForientation = 0
                elif myframe == 1:
                    ORForientation = 2
                elif myframe == 2:
                    ORForientation = 1
                else:
                    continue
            else:
                continue
        # Determine nucleotide window around which to get sequence
        wholeseqstart = 0
        wholeseqend = 0
        if strand == ('+'):
            if (peptideends - (exonends-(nAAs*3-3)+ORForientation)) >= (nAAs*3):
                wholeseqstart = exonends-(nAAs*3-3)+ORForientation
                wholeseqend = peptideends
            else:
                continue
        else:
            if ((exonends-(nAAs*3-3)-ORForientation)-peptideends) >= (nAAs*3):
                wholeseqstart = peptideends
                wholeseqend = exonstarts+(nAAs*3-3)-ORForientation
            else:
                continue
        #    if strand == ('+'):
        #        wholeseqstart = exonstarts - exonframes - (nAAs*3)
        #        wholeseqend = peptideends
        #    else:
        #        wholeseqstart = peptideends 
        #        wholeseqend =  exonstarts + exonframes + (nAAs*3)
                # Get genomic sequence
        if strand == ('+'):
            loc = chrom + ':' + str(wholeseqstart - 1) + '-' + str(wholeseqend)
        else:
            loc = chrom + ':' + str(wholeseqstart - 1) + '-' + str(wholeseqend)
        loc = '-seq='+loc
                #print loc
        twobitoutput = (subprocess.check_output(['/data/11000039/e0149673/scratch/bin/anaconda2/bin/twoBitToFa',loc,'/data/11000039/e0149673/scratch/bin/library/hg19.2bit','stdout']))
                # Parse output
        seqlist = twobitoutput.split('\n')
        headerline = seqlist[0]+"|"+tpmval
        seqlist = seqlist[1:len(seqlist)-1]
        sequence = ''.join(str(elem) for elem in seqlist)
        sequence = sequence.upper()
        #print sequence
        # Reverse complement if it's on the negative strand
        if strand == ('-'):
            sequence = str(Seq(sequence).reverse_complement())
        # Manually translate sequence
        peptide = manualTranslate(sequence)
        # Check to make sure peptide is at least "length" AAs long, and if so, write to output file
        if len(peptide) < nAAs:
            continue
        else:
            loc = chrom + ':' + str(wholeseqstart) + ';' + str(wholeseqstart+(len(peptide)*3)) + ';' + str(wholeseqend)
            if strand == ('-'):
                loc = chrom + ':' + str(wholeseqend) + ';' + str(wholeseqend-(len(peptide)*3)) + ';' + str(wholeseqstart)
            #loc = '-seq='+loc
            newheaderline = '>s='
            outfile.write(newheaderline+loc+','+strand+','+tpmval+','+gene_name)
            #outfile.write(newheaderline+loc+','+strand+','+tpmval+','+gene_name+','+file_name)
            outfile.write(peptide+'\n')
            headermapfile.write(newheaderline+'\t'+headerline+'\n')
            #print peptide
    return
    
# Main function that processes command line input and calls other functions
def main():
    # Check to make sure we have the right number of inputs
    if len(sys.argv) != 4:
        print('Error: incorrect number of inputs.')
        print('Please input the AA window you want, an outfile path, and a polyA site file.')
        sys.exit()
    # Store inputs
#        kmafile = sys.argv[1]
    window = int(sys.argv[1])
    outpath = sys.argv[2]
#    uniqueIntrons = open("/data2/renxi/bin/retained-intron-neoantigen-pipeline/text.txt", "r")
    #_, indices = np.unique(chromlocs, return_index=True)
    #uniquelocs = chromlocs[indices,:]
#        uniqueIntrons = np.core.defchararray.strip(uniquelocs, '"')  # Strip "s from locations
        # Create unique intron output file
        #uniqueIntrons,TPMvals = createUniqueIntronList(kmafile, outpath)
#        window = '9'
        # Create nucleotide sequences file
    getSeqs(window, outpath)

if __name__ == '__main__':
    main()
