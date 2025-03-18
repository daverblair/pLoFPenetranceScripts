import pandas as pd
import numpy as np
import pysam
from cyvcf2 import VCF
from Bio.Seq import Seq
from Bio.Data.CodonTable import TranslationError



def NormalizeVariantCsq(csq_list):
	"""Based on of annotations, returns one of the three pLOF classes: FRAMESHIFT, STOPGAINED, and SPLICE_CHANGE. 
	
	Args:
	    csq_list (list): list of annotations for variant
	
	Returns:
	    Str: pLOF class
	
	Raises:
	    ValueError: list of annotations does not fit into a pLOF class
	"""
	if 'frameshift_variant' in csq_list:
		return 'FRAMESHIFT'
	elif 'stop_gained' in csq_list:
		return 'STOP_GAINED'
	elif 'splice' in set([x.split('_')[0] for x in csq_list]):
		return 'SPLICE_CHANGE'
	else:
		print(csq_list)
		raise ValueError('Variant does not fit into a pLOF class')


def AlternateVariantIndex(v_id):
	"""Normalizes variant in the following string format: CHROM_POS_REF_ALT
	
	Args:
	    v_id (str): CHROM_POS_REF_ALT
	
	Returns:
	    Str: normalized variant
	"""
	chrom,pos,ref,alt=v_id.split('_')
	if (len(ref)>1) and (len(ref)==len(alt)):
		num_substitutions=sum(1 for a, b in zip(ref, alt) if a != b)
		substitution_positions=[i for i in range(len(ref)) if ref[i] != alt[i]]
		if num_substitutions==1:
			#all nucleotides but 1 are unchanged, so drop them and adjust the position
			substitution_position=substitution_positions[0]
			new_variant_pos=int(pos)+substitution_position
			new_variant_ref=ref[substitution_position]
			new_variant_alt=alt[substitution_position]
		else:
			#more than 1 nucleotide differs, so keep all the different ones
			new_variant_pos = int(pos)+min(substitution_positions)
			new_variant_ref=ref[min(substitution_positions):max(substitution_positions)]
			new_variant_alt=alt[min(substitution_positions):max(substitution_positions)]

		alternate_id='_'.join([chrom,str(new_variant_pos),new_variant_ref,new_variant_alt])

		return alternate_id,{'CHROM':chrom,'POS':new_variant_pos,'REF':new_variant_ref,'ALT':new_variant_alt}
	elif (len(ref)>1) and (len(ref)<len(alt)):
		#insertion
		for i,nuc in enumerate(ref):
			if nuc!=alt[i]:
				break
		new_variant_pos=int(pos)+i
		new_variant_ref=ref[i]
		new_variant_alt=alt[i:i+(len(alt)-len(ref))+1]
		alternate_id='_'.join([chrom,str(new_variant_pos),new_variant_ref,new_variant_alt])
		return alternate_id,{'CHROM':chrom,'POS':new_variant_pos,'REF':new_variant_ref,'ALT':new_variant_alt}

	elif (len(ref)>1) and (len(ref)>len(alt)):
		#deletion
		for i,nuc in enumerate(alt):
			if nuc!=ref[i]:
				break
		new_variant_pos=int(pos)+i

		new_variant_ref=ref[i:i+(len(ref)-len(alt))+1]
		new_variant_alt=alt[i]
		alternate_id='_'.join([chrom,str(new_variant_pos),new_variant_ref,new_variant_alt])
		return alternate_id,{'CHROM':chrom,'POS':new_variant_pos,'REF':new_variant_ref,'ALT':new_variant_alt}
	else:
		return v_id,{'CHROM':chrom,'POS':int(pos),'REF':ref,'ALT':alt}




def ComputeCDSLength(exons,utr_5prime,utr_3prime):
	"""Given a set of exons and the 5' and 3' UTRs, return the coding sequence length
	
	Args:
	    exons (list of tuples): list of exons
	    utr_5prime (list of tuples): list of 5' utr exon tuples
	    utr_3prime (list of tuples): list of 3' utr exon tuples
	
	Returns:
	    int: coding nucleotide length ()
	"""
	return sum([(x[1]-x[0])+1 for x in exons]) - sum([(x[1]-x[0])+1 for x in utr_5prime]) - sum([(x[1]-x[0])+1 for x in utr_3prime])

def ComputeAALength(cds_length):
	"""Returns the AA length from the coding sequence
	
	Args:
	    cds_length (int): Lenfth of conding sequence
	
	Returns:
	    int: length of AA sequence
	"""
	assert cds_length % 3 == 0, "CDS is not a multiple of 3"
	return (cds_length / 3)-1 # remove 1 for the stop codon

def BuildSequences(transcript,tx_start,tx_end,exon_list,utr_5prime,utr_3prime,strand):
	"""Attempts to build both the cDNA and AA sequence from the full transcript sequence, the exons, and the 5'/3' UTRs. Note, many non-canonical transcripts have non-sensical AA translations. This function tries to recover a coherent translation for downstream analyses of sequence context, which may differ from what's in the reference databases. Warnings are produced for various errors, but in general, errors are forced through the function. 
	
	Args:
	    transcript (str): string of transcript nucleotides
	    tx_start (int): genome position of transcript start
	    tx_end (int): genome position of transcript start
	    exon_list (list of tuples): List of exons in transcript
	    utr_5prime (list of tuples): List of 5' UTR regions
	    utr_3prime (list of tuples): List of 3' UTR regions
	    strand (int): strand, must be 1 or -1
	
	Returns:
	    (str,str): Returns the cDNA and AA sequences. Note, for canonical/MANE select transcripts, these will agree with database versions. For non-canonical transcripts of varying quality, there will be discrepancies. 
	"""
	assert len(transcript)==(tx_end-tx_start)+1,"Transcript sequence length does not match the tx_start and tx_end."
	if strand==1:
		transformed_exons=[((x[0]-tx_start)+1,(x[1]-tx_start)+1)for x in exon_list]
	else:
		transformed_exons=[(len(transcript)-(x[1]-tx_start),len(transcript)-(x[0]-tx_start)) for x in exon_list]
	#cut out introns

	utr_5prime_len=sum([(x[1]-x[0])+1 for x in utr_5prime])
	utr_3prime_len=sum([(x[1]-x[0])+1 for x in utr_3prime])
	cDNA_untrimmed = Seq('').join([transcript[x[0]-1:x[1]] for x in transformed_exons])
	cDNA = cDNA_untrimmed[utr_5prime_len:len(cDNA_untrimmed)-utr_3prime_len]
	try:
		AA=cDNA.translate(cds=True)
		if '*' in AA:
			# if a stop codon is in the AA sequence, then return the truncated AA sequence. Note, it's possible that this stop codon is read through, but difficult to be certain. 
			AA=AA[:AA.find('*')]
	# if cDNA doesn't translate according to BioSeq's expectatons...
	except TranslationError:
		#translate in error mode
		error_AA=str(cDNA.translate())
		start_codon=cDNA[0:3]
		#non-canonical start codon, use it
		if start_codon in set(['ACG','CTG', 'GTG', 'TTG']):
			print('Warning: Transcript with start site {0:d} has non-canonical start codon.'.format(tx_start)) 
			error_AA=list(error_AA)
			error_AA[0]='M'
			error_AA=''.join(error_AA)
		else:
			# no apparent non-canonical start, so will just use the first methionine, assuming that there is some 5' UTR error
			if 'M' in error_AA:
				error_AA=error_AA[error_AA.find('M'):]
				print('Warning: transcript with start site {0:d} does not have a start codon a position 1. Please check.'.format(tx_start))
			else:
				# at this point, very unlikely that this represents a PC coding transcript
				print('Warning: transcript with start site {0:d} has no start codon. Please check.'.format(tx_start))
		
		if error_AA[-1]!='*':
			# should end with a stop codon, but does not
			print('Transcript with start site {0:d} has no terminal stop codon. Please check.'.format(tx_start))
			# assuming that the 3' UTR is incorrect, will truncate at the first stop codon
			if error_AA.count('*')>0:
				error_AA=error_AA[:error_AA.find('*')]
			else:
				# no stop in transcript. Will return, but unlikely that this is a PC transcript. 
				print('Transcript with start site {0:d} has no stop codon. Please check.'.format(tx_start))
		else:
			# ends with a stop codon, which is good
			error_AA=error_AA[:-1]
			#checking for other stop codons mid-transcript
			if error_AA.count('*')>0:
				error_AA=error_AA[:error_AA.find('*')]
				print('Transcript with start site {0:d} has extra stop codon. Please check.'.format(tx_start))
		AA=Seq(error_AA)
	return cDNA,AA


def ReturnFirstCodingExon(exon_list,utr_5prime,utr_3prime,strand):
	"""Returns the first coding exon in a transcript. 
	
	Args:
	    exon_list (list of tuples): List of exons in transcript
	    utr_5prime (list of tuples): List of 5' UTR regions
	    utr_3prime (list of tuples): List of 3' UTR regions
	    strand (int): strand, must be 1 or -1
	
	Returns:
	    (tuple, tuple): first coding exon, coding portion of exon
	

	"""
	first_found=False
	utr_5prime_overlaps=False
	if len(utr_5prime)>0:
		for first_coding_exon in exon_list:
			if (strand==1):
				if(utr_5prime[-1][0]>=first_coding_exon[0]) and (utr_5prime[-1][1]<first_coding_exon[1]):
					first_found=True
					utr_5prime_overlaps=True
					break
			else:
				if (utr_5prime[-1][1]<=first_coding_exon[1]) and (utr_5prime[-1][0]>first_coding_exon[0]):
					first_found=True
					utr_5prime_overlaps=True
					break
	else:
		first_coding_exon=exon_list[0]
		first_found=True


	if first_found==False:
		#check if one perfectly overlaps
		for first_coding_exon in exon_list:
			if (strand==1):
				if (utr_5prime[-1][0]>=first_coding_exon[0]) and (utr_5prime[-1][1]<=first_coding_exon[1]):
					first_found=True
					break
			else:
				if (utr_5prime[-1][1]<=first_coding_exon[1]) and (utr_5prime[-1][0]>=first_coding_exon[0]):
					first_found=True
					break
		first_coding_exon=exon_list[exon_list.index(first_coding_exon)+1]

	if first_found==False:
		raise ValueError('Unable to find first exon')

	if utr_5prime_overlaps==True:
		if (strand==1):
			coding_portion=[utr_5prime[-1][1]+1,first_coding_exon[1]]
		else:
			coding_portion=[first_coding_exon[0],utr_5prime[-1][0]-1]
	else:
		coding_portion=[first_coding_exon[0],first_coding_exon[1]]

	#now trim off 3'UTR if single exon gene
	utr_3_prime_overlaps=False
	if len(utr_3prime)>0:
		if strand==1:
			if (utr_3prime[0][0]>first_coding_exon[0]) and (utr_3prime[0][1]<=first_coding_exon[1]):
				utr_3_prime_overlaps=True
		else:
			if (utr_3prime[0][1]<first_coding_exon[1]) and (utr_3prime[0][0]>=first_coding_exon[0]):
				utr_3_prime_overlaps=True

	if utr_3_prime_overlaps==True:
		if (strand==1):
			coding_portion=[coding_portion[0],utr_3prime[0][0]-1]
		else:
			coding_portion=[utr_3prime[0][1]+1,coding_portion[1]]

	return first_coding_exon,tuple(coding_portion)


def ReturnLastCodingExon(exon_list,utr_5prime,utr_3prime,strand):
	"""Returns the last coding exon in a transcript. 
	
	Args:
	    exon_list (list of tuples): List of exons in transcript
	    utr_5prime (list of tuples): List of 5' UTR regions
	    utr_3prime (list of tuples): List of 3' UTR regions
	    strand (int): strand, must be 1 or -1
	
	Returns:
	    (tuple, tuple): last coding exon, coding portion of exon
	
	"""

	last_found=False
	utr_3_prime_overlaps=False
	if len(utr_3prime)>0:
		for last_coding_exon in exon_list:
			if (strand==1):
				if (utr_3prime[0][0]>last_coding_exon[0]) and (utr_3prime[0][1]<=last_coding_exon[1]):
					last_found=True
					utr_3_prime_overlaps=True
					break
			else:
				if (utr_3prime[0][1]<last_coding_exon[1]) and (utr_3prime[0][0]>=last_coding_exon[0]):
					last_found=True
					utr_3_prime_overlaps=True
					break
		if last_found==False:
			for last_coding_exon in exon_list:
				if (strand==1):
					if (utr_3prime[0][0]>=last_coding_exon[0]) and (utr_3prime[0][1]<=last_coding_exon[1]):
						last_found=True
						break
				else:
					if (utr_3prime[0][1]<=last_coding_exon[1]) and (utr_3prime[0][0]>=last_coding_exon[0]):
						last_found=True
						break
			last_coding_exon=exon_list[exon_list.index(last_coding_exon)-1]
	else:
		last_coding_exon=exon_list[-1]
		last_found=True

	if last_found==False:
		raise ValueError('Unable to find last exon')

	if utr_3_prime_overlaps==True:
		if (strand==1):
			coding_portion=[last_coding_exon[0],utr_3prime[0][0]-1]
		else:
			coding_portion=[utr_3prime[0][1]+1,last_coding_exon[1]]
	else:
		coding_portion=[last_coding_exon[0],last_coding_exon[1]]


	#now trim off 5'UTR if single exon gene
	utr_5prime_overlaps=False
	if len(utr_5prime)>0:
		if strand==1:
			if  (utr_5prime[-1][0]>=last_coding_exon[0]) and (utr_5prime[-1][1]<last_coding_exon[1]):
				utr_5prime_overlaps=True
		else:
			if (utr_5prime[-1][1]<=last_coding_exon[1]) and (utr_5prime[-1][0]>last_coding_exon[0]):
				utr_5prime_overlaps=True

	if utr_5prime_overlaps==True:
		if (strand==1):
			coding_portion=[utr_5prime[-1][1]+1,coding_portion[1]]
		else:
			coding_portion=[coding_portion[0],utr_5prime[-1][0]-1]

	return last_coding_exon,coding_portion
			


def ReturnLastExonAARange(AA_length,exon_list,utr_5prime,utr_3prime,strand):
	"""Returns the AA range for the last coding exon
	
	Args:
	    AA_length (int): length of AA sequence
	    exon_list (list of tuples): List of exons in transcript
	    utr_5prime (list of tuples): List of 5' UTR regions
	    utr_3prime (list of tuples): List of 3' UTR regions
	    strand (int): strand, must be 1 or -1
	
	Returns:
	    (tuple): last coding exon AA range
	"""
	last_coding_exon,last_coding_portion=ReturnLastCodingExon(exon_list,utr_5prime,utr_3prime,strand)
	#subtract off the stop codon
	num_aa=np.ceil((last_coding_portion[1]-last_coding_portion[0]+1)/3)-1
	return (int(AA_length-num_aa+1),AA_length)

def ReturnFirstExonAARange(AA_length,exon_list,utr_5prime,utr_3prime,strand):
	"""Returns the AA range for the first coding exon
	
	Args:
	    AA_length (int): length of AA sequence
	    exon_list (list of tuples): List of exons in transcript
	    utr_5prime (list of tuples): List of 5' UTR regions
	    utr_3prime (list of tuples): List of 3' UTR regions
	    strand (int): strand, must be 1 or -1
	
	Returns:
	    (tuple): first coding exon AA range
	"""

	first_coding_exon,first_coding_portion=ReturnFirstCodingExon(exon_list,utr_5prime,utr_3prime,strand)
	last_coding_exon,last_coding_portion=ReturnLastCodingExon(exon_list,utr_5prime,utr_3prime,strand)
	num_aa=np.ceil((first_coding_portion[1]-first_coding_portion[0]+1)/3)
	#subtract off the stop codon if the first is the same as the last
	if last_coding_exon==first_coding_exon:
		num_aa-=1
	return (1,int(num_aa))


def ParseSpliceAI(file_handle):
	"""Parses a splice AI VCF file.
	
	Args:
	    file_handle (str): filename
	
	Returns:
	    pd.DataFrame: Dataframe containing the splice AI results. Each variant-transcript pair has a unique line, indexed by variants (ie redundant index)
	"""
	spliceai_vcf = VCF(file_handle)
	header_id_dict={}
	for header_entry in spliceai_vcf.header_iter():
		if header_entry.type=='INFO':
			header_id_dict[header_entry.info()['ID']]=header_entry.info()['Description']

	SpliceAI_Columns=header_id_dict['SpliceAI'].strip('"').split('Format:')[1].strip().split('|')
	output_data={'VARIANT_ID':[],'TX_ID':[],'SPLICE_AI':[]}

	for variant in spliceai_vcf:
		try:
			info=variant.INFO.get('SpliceAI').split(',')
			num_transcripts=int(len(info)/len(variant.ALT))	
			for i,alt in enumerate(variant.ALT):
				for j in range(num_transcripts):
					if variant.CHROM[0:3]!='chr':
						v_id='chr'+variant.CHROM+'_'+str(variant.POS)+'_'+variant.REF+'_'+alt
					else:
						v_id=variant.CHROM+'_'+str(variant.POS)+'_'+variant.REF+'_'+alt
					symbol=info[i*num_transcripts+j].split('|')[1]
					parsed_info=dict(zip(SpliceAI_Columns[2:],info[i*num_transcripts+j].split('|')[2:]))
					for key,value in parsed_info.items():
						if key.split('_')[0]=='DS':
							try:
								parsed_info[key]=float(value)
							except ValueError:
								parsed_info[key]=pd.NA
						else:
							try:
								parsed_info[key]=int(value)
							except ValueError:
								parsed_info[key]=pd.NA

					output_data['VARIANT_ID']+=[v_id]
					output_data['SPLICE_AI']+=[parsed_info]
					output_data['TX_ID']+=[symbol]
		except AttributeError:
			pass
	output_data=pd.DataFrame(output_data)
	output_data.set_index('VARIANT_ID',inplace=True)
	return output_data


def ReturnFirstImpactedExon(splice_change_type,genome_position,strand,exons):
	"""Returns the first impacted exon of a splice change. For a donor change, it's the exon downstream. For an acceptor change, it's the same exon with the acceptor. The impacted exon is chosen based on the genomic position of the change and it's type. 
	
	Args:
	    splice_change_type (str): type of splice change
	    genome_position (int): genomic position of the change
	    strand (int): 1 or -1
	    exons (list of tuples): list of exons
	
	Returns:
	    tuple: impacted exon
	"""
	if splice_change_type.split('_')[0]=='DONOR':
		#targeted exon is the one with the donor site
		if strand==-1:
			all_possible_donor_sites=[x[0] for x in exons]
		else:
			all_possible_donor_sites=[x[1] for x in exons]

		genome_distance=[np.abs(genome_position-x) for x in all_possible_donor_sites]
		#donor site mutation impacts the downstream exon
		impacted_exon=exons[genome_distance.index(min(genome_distance))+1]
	else:
		if strand==-1:
			all_possible_acceptor_sites=[x[1] for x in exons]
		else:
			all_possible_acceptor_sites=[x[0] for x in exons]
		genome_distance=[np.abs(genome_position-x) for x in all_possible_acceptor_sites]
		#acceptor site mutation impacts the closest exon
		impacted_exon=exons[genome_distance.index(min(genome_distance))]
	return impacted_exon



def ReturnPrimarySpliceInfo(splice_ai_data,variant_pos,strand,exons,utr_5prime,utr_3prime):
	"""Given a dictionary of SpliceAI results, the position of the variant, and transcript information, the following information is returned:

	1) the highest scoring splice change type (typically donor or acceptor loss but sometimes a gain is highest)
	2) the score associated with change
	3) If translated, the fraction of the coding nucleotides predicted to be impacted by the change

	Args:
	    splice_ai_data (dict): Dictionary with the following entries: {'DS_AG','DS_AL','DS_DG','DS_DL','DP_AG','DP_AL','DP_DG':,'DP_DL'}
	    variant_pos (int): genomic position
	    strand (int): 1 or -1
	    exons (list of tuples): transcript exons
	    utr_5prime (list of tuples): tuples defining 5' UTR
	    utr_3prime (list of tuples): tuples defining 3' UTR
	
	Returns:
	    : See above
	"""
	conv_key={'DS_AG':'ACC_GAIN','DS_AL':'ACC_LOSS','DS_DG':'DONOR_GAIN','DS_DL':'DONOR_LOSS','DP_AG':'ACC_GAIN','DP_AL':'ACC_LOSS','DP_DG':'DONOR_GAIN','DP_DL':'DONOR_LOSS'}
	score_table=pd.DataFrame([],columns=['SCORE','GENOME_POS'],index=['ACC_GAIN','ACC_LOSS','DONOR_GAIN','DONOR_LOSS'])
	for key,value in splice_ai_data.items():
		column='SCORE' if key.split('_')[0]=='DS' else 'GENOME_POS'
		if column=='GENOME_POS':
			try:
				genome_pos=variant_pos+value
			except ValueError:
				genome_pos=pd.NA

			score_table.loc[conv_key[key],column]=genome_pos
		else:
			score_table.loc[conv_key[key],column]=value

	# check if in-frame with donor loss


	if score_table['SCORE'].sum()>0.0:
		max_splice_score=score_table.SCORE.max()
		splice_type=score_table.SCORE.idxmax()
		impacted_position=score_table.loc[score_table.SCORE.idxmax()]['GENOME_POS']

		#now, figure out the fraction impacted
		first_impacted_exon=ReturnFirstImpactedExon(splice_type,impacted_position,strand,exons)
		first_impacted_exon_number=exons.index(first_impacted_exon)+1
		affected_plus_downstream_exons=exons[first_impacted_exon_number-1:]

		first_coding_exon,first_coding_portion=ReturnFirstCodingExon(exons,utr_5prime,utr_3prime,strand)
		last_coding_exon,last_coding_portion=ReturnLastCodingExon(exons,utr_5prime,utr_3prime,strand)

		if first_coding_exon in affected_plus_downstream_exons:
			first_exon_index = affected_plus_downstream_exons.index(first_coding_exon)
			affected_plus_downstream_exons=affected_plus_downstream_exons[first_exon_index:]
			affected_plus_downstream_exons[0]=tuple(first_coding_portion)
		if last_coding_exon in affected_plus_downstream_exons:
			last_exon_index = affected_plus_downstream_exons.index(last_coding_exon)
			affected_plus_downstream_exons=affected_plus_downstream_exons[:last_exon_index+1]
			affected_plus_downstream_exons[-1]=tuple(last_coding_portion)

		if (strand==1) and ((affected_plus_downstream_exons[-1][1]<first_coding_portion[0]) or (affected_plus_downstream_exons[0][0]>last_coding_portion[1])):
			frac_coding_nuc_impacted=0.0
		elif (strand==-1) and ((affected_plus_downstream_exons[-1][0]>first_coding_portion[1]) or (affected_plus_downstream_exons[0][1]<last_coding_portion[0])):
			frac_coding_nuc_impacted=0.0
		else:

			num_coding_nucleotides_lost=sum([(x[1]-x[0])+1 for x in affected_plus_downstream_exons]) 
			#subtract 1 to account for stop codon at the end of the nucleotide sequence
			frac_coding_nuc_impacted=num_coding_nucleotides_lost/ComputeCDSLength(exons,utr_5prime,utr_3prime)

	else:
		if pd.isna(score_table['SCORE']).sum()==4:
			max_splice_score=0.0
			splice_type='INDETERMINATE'
			frac_coding_nuc_impacted=0.0
			
		else:
			max_splice_score=0.0
			splice_type='INDETERMINATE'
			frac_coding_nuc_impacted=0.0
	return splice_type,max_splice_score,frac_coding_nuc_impacted



def PredictAcceptorLossRescue(splice_ai_data,variant_pos,strand,exons,utr_5prime,utr_3prime,aa_sequence,met_rescue_threshold=0.9):
	"""Given an acceptor loss, predicts possible rescue events. These include:

	1) Original splice site intact: the change simply rescues/impacts the same normal splice site without any clear splicing effect. 

	2) Non-coding exon: if the acceptor impacts a non-coding exon, then it's unlikely to have an effect. 

	3) Last exon: an acceptor loss in the last exon is more likely to be tolerated. 
	
	4) Cryptic Rescue: an inframe gain rescues the loss

	5) Exon skipping: An inframe exon is skipped

	6) Intron retained: an inframe intron can be retained

	7) Methionine resuce: The first coding exon is impacted by the loss, but there is a methionine downstream within the first 1.0 - met_rescue_threshold fraction of the protein. 

	
	Args:
	    splice_ai_data (dict): Dictionary with the following entries: {'DS_AG','DS_AL','DS_DG','DS_DL','DP_AG','DP_AL','DP_DG':,'DP_DL'}
	    variant_pos (int): genomic position
	    strand (int): 1 or -1
	    exons (list of tuples): transcript exons
	    utr_5prime (list of tuples): tuples defining 5' UTR
	    utr_3prime (list of tuples): tuples defining 3' UTR
	    aa_sequence (str): AA sequence
	    met_rescue_threshold (float, optional): threshold for methionine rescue. If the impacted exon is skipped and is the first exon, opportunity for methionine rescue. At least met_rescue_threshold fraction of the AA sequence must remain. 
	
	Returns:
	    list of tuples: List of tuples for all possible rescues. Each tuple contains: (label, confidence score, fraction coding nucs impacted)
	"""
	conv_key={'DS_AG':'ACC_GAIN','DS_AL':'ACC_LOSS','DS_DG':'DONOR_GAIN','DS_DL':'DONOR_LOSS','DP_AG':'ACC_GAIN','DP_AL':'ACC_LOSS','DP_DG':'DONOR_GAIN','DP_DL':'DONOR_LOSS'}
	score_table=pd.DataFrame([],columns=['SCORE','GENOME_POS'],index=['ACC_GAIN','ACC_LOSS','DONOR_GAIN','DONOR_LOSS'])
	for key,value in splice_ai_data.items():
		column='SCORE' if key.split('_')[0]=='DS' else 'GENOME_POS'
		if column=='GENOME_POS':
			try:
				genome_pos=variant_pos+value
			except ValueError:
				genome_pos=pd.NA

			score_table.loc[conv_key[key],column]=genome_pos
		else:
			score_table.loc[conv_key[key],column]=value

	acc_loss_pos = score_table.loc['ACC_LOSS','GENOME_POS']
	acc_gain_pos = score_table.loc['ACC_GAIN','GENOME_POS']
	target_exon=ReturnFirstImpactedExon('ACC_LOSS',acc_loss_pos,strand,exons)
	upstream_exon = exons[exons.index(target_exon)-1]
	first_coding_exon,first_coding_portion=ReturnFirstCodingExon(exons,utr_5prime,utr_3prime,strand)
	last_coding_exon,last_coding_portion=ReturnLastCodingExon(exons,utr_5prime,utr_3prime,strand)
	cds_length=ComputeCDSLength(exons,utr_5prime,utr_3prime)

	if (strand==1):
		prior_acceptor_site = target_exon[0]
	else:
		prior_acceptor_site=target_exon[1]

	all_flags = []

	# Original splice site intact
	if prior_acceptor_site!=acc_loss_pos:
		#the original site is intact, so the score is the difference between the site (presumably 1) and the score of the gain. Loss at any other site is irrelevant.
		all_flags+=[('ORIG_SPLICESITE_INTACT',1.0-score_table.loc['ACC_GAIN','SCORE'],0.0)]

	elif score_table.loc['ACC_GAIN','GENOME_POS']==prior_acceptor_site:
		#the original site is the gain, so return the score. 
		all_flags+=[('ORIG_SPLICESITE_INTACT',score_table.loc['ACC_GAIN','SCORE'],0.0)]

	#outside coding region
	if (strand==1):
		if (prior_acceptor_site>last_coding_exon[1]) or (prior_acceptor_site<first_coding_exon[0]):
			all_flags+=[('OUTSIDE_CODING_REGION',1.0,0.0)]
	else:
		if (prior_acceptor_site<last_coding_exon[0]) or (prior_acceptor_site>first_coding_exon[1]):
			all_flags+=[('OUTSIDE_CODING_REGION',1.0,0.0)]

	#last exon
	if (last_coding_exon==target_exon) and (first_coding_exon!=last_coding_exon):
		last_coding_frac_lost = (last_coding_portion[1]-last_coding_portion[0]+1)/cds_length
		all_flags+=[('LAST_CODING_EXON',1.0,last_coding_frac_lost)]

	# Cryptic Rescue (with respect to original transcript site)
	nucleotide_difference = np.abs(acc_gain_pos-prior_acceptor_site)

	# inframe cryptic gain
	if ((nucleotide_difference % 3)==0) and (score_table.loc['ACC_GAIN','SCORE']>0.0):
		cryptic_rescue_score  = score_table.loc['ACC_GAIN','SCORE']
		if (strand==1):
			if (acc_gain_pos > prior_acceptor_site):
				cryptic_rescue_frac_coding_nuc_impacted= nucleotide_difference/cds_length
			else:
				cryptic_rescue_frac_coding_nuc_impacted=nucleotide_difference/(cds_length+nucleotide_difference)
		else:
			if (acc_gain_pos < prior_acceptor_site):
				cryptic_rescue_frac_coding_nuc_impacted= nucleotide_difference/cds_length
			else:
				cryptic_rescue_frac_coding_nuc_impacted=nucleotide_difference/(cds_length+nucleotide_difference)
		all_flags+=[('CRYPTIC_RESCUE',cryptic_rescue_score,cryptic_rescue_frac_coding_nuc_impacted)]


	#internal exon skipping
	if (target_exon!=first_coding_exon) and (target_exon!=last_coding_exon):
		target_exon_size = (target_exon[1]-target_exon[0])+1
	elif (target_exon==first_coding_exon):
		target_exon_size = (first_coding_portion[1]-first_coding_portion[0]+1)
	else:
		target_exon_size = (last_coding_portion[1]-last_coding_portion[0]+1)

	if (target_exon_size % 3)==0:
		exon_skip_score = 1.0-score_table.loc['ACC_GAIN','SCORE']
		exon_skip_frac_coding_nuc_impacted = target_exon_size/cds_length
		all_flags+=[('INFRAME_EXON_SKIP',exon_skip_score,exon_skip_frac_coding_nuc_impacted)]


	#Here, we have to lose the upstream donor for intron retention. This is only possible assuming the acceptor isn't in the first exon. 
	if (target_exon!=exons[0]):
		if (strand==1):
			intron_size =target_exon[0]-upstream_exon[1]-1
		else:
			intron_size =upstream_exon[0]-target_exon[1]-1

		if ((intron_size % 3)==0):
			if (strand==1):
				upstream_donor_position=upstream_exon[1]
			else:
				upstream_donor_position=upstream_exon[0]

			if score_table.loc['DONOR_LOSS','GENOME_POS']==upstream_donor_position:
				intron_retention_score = score_table.loc['DONOR_LOSS','SCORE']
				intron_retention_frac_coding_nuc_impacted = intron_size/(ComputeCDSLength(exons,utr_5prime,utr_3prime)+intron_size)
			else:
				intron_retention_score=0.0
				intron_retention_frac_coding_nuc_impacted=intron_size/(ComputeCDSLength(exons,utr_5prime,utr_3prime)+intron_size)
			all_flags+=[('INFRAME_INTRON_RETENTION',intron_retention_score,intron_retention_frac_coding_nuc_impacted)]

	#methionine rescue of first exon skip
	if (target_exon==first_coding_exon):
		first_exon_range = ReturnFirstExonAARange(len(aa_sequence),exons,utr_5prime,utr_3prime,strand)
		if len(aa_sequence[first_exon_range[1]:])/len(aa_sequence) >= met_rescue_threshold:
			#note, off by one indexing is intentional
			new_sequence = aa_sequence[first_exon_range[1]:]
			first_methioninine = new_sequence.find('M')
			new_sequence=new_sequence[new_sequence.find('M'):]
			if len(new_sequence)/len(aa_sequence) >= met_rescue_threshold:
				possible_met_rescue=True
				met_rescue_nuc_impacted = 1.0-(len(aa_sequence[first_exon_range[1]:])/len(aa_sequence))
				all_flags+=[('POSSIBLE_MET_RESCUE',1.0,met_rescue_nuc_impacted)]
	return all_flags



def PredictAcceptorGainRescue(splice_ai_data,variant_pos,strand,exons,utr_5prime,utr_3prime):
	"""
	Given an acceptor gain, predict rescues as follows:

	1) Original splice site intact: it simply rescues/impacts the same normal splice site

	2) Non-coding exon: if the acceptor impacts a non-coding exon, then it's unlikely to have an effect. 

	3) Last exon: an aberrant acceptor gain in the last exon is more likely to be tolerated.
	
	4) Cryptic Rescue: the inframe gain rescues a loss at the normal splice site

	Note, because an acceptor gain induces a new splice site, it cannot be overcome by exon skipping, intron retention, or methionine rescue (assuming that it occurs in isolation). Theoretically, if accompanied by an acceptor loss, skipping could occur. But since the gain is more confident than the loss, will defer this possibility for now. 

	
	Args:
	    splice_ai_data (dict): Dictionary with the following entries: {'DS_AG','DS_AL','DS_DG','DS_DL','DP_AG','DP_AL','DP_DG':,'DP_DL'}
	    variant_pos (int): genomic position
	    strand (int): 1 or -1
	    exons (list of tuples): transcript exons
	    utr_5prime (list of tuples): tuples defining 5' UTR
	    utr_3prime (list of tuples): tuples defining 3' UTR
	
	Returns:
	     list of tuples: List of tuples for all possible rescues. Each tuple contains: (label, confidence score, fraction coding nucs impacted)
	"""
	conv_key={'DS_AG':'ACC_GAIN','DS_AL':'ACC_LOSS','DS_DG':'DONOR_GAIN','DS_DL':'DONOR_LOSS','DP_AG':'ACC_GAIN','DP_AL':'ACC_LOSS','DP_DG':'DONOR_GAIN','DP_DL':'DONOR_LOSS'}
	score_table=pd.DataFrame([],columns=['SCORE','GENOME_POS'],index=['ACC_GAIN','ACC_LOSS','DONOR_GAIN','DONOR_LOSS'])
	for key,value in splice_ai_data.items():
		column='SCORE' if key.split('_')[0]=='DS' else 'GENOME_POS'
		if column=='GENOME_POS':
			try:
				genome_pos=variant_pos+value
			except ValueError:
				genome_pos=pd.NA

			score_table.loc[conv_key[key],column]=genome_pos
		else:
			score_table.loc[conv_key[key],column]=value

	acc_loss_pos = score_table.loc['ACC_LOSS','GENOME_POS']
	acc_gain_pos = score_table.loc['ACC_GAIN','GENOME_POS']
	target_exon=ReturnFirstImpactedExon('ACC_GAIN',acc_gain_pos,strand,exons)
	upstream_exon = exons[exons.index(target_exon)-1]
	first_coding_exon,first_coding_portion=ReturnFirstCodingExon(exons,utr_5prime,utr_3prime,strand)
	last_coding_exon,last_coding_portion=ReturnLastCodingExon(exons,utr_5prime,utr_3prime,strand)
	cds_length=ComputeCDSLength(exons,utr_5prime,utr_3prime)

	if (strand==1):
		prior_acceptor_site = target_exon[0]
	else:
		prior_acceptor_site=target_exon[1]

	all_flags = []

	# Original splice site intact
	if acc_gain_pos==prior_acceptor_site:
		#the original site is the gain, so return the score. 
		all_flags+=[('ORIG_SPLICESITE_INTACT',score_table.loc['ACC_GAIN','SCORE'],0.0)]

	elif prior_acceptor_site!=acc_loss_pos:
		#the original site is intact, so the score is the difference between the site (presumably 1) and the score of the gain. Loss at any other site is irrelevant.
		all_flags+=[('ORIG_SPLICESITE_INTACT',1.0-score_table.loc['ACC_GAIN','SCORE'],0.0)]

	#gain occurs outside coding region
	if (strand==1):
		if (acc_gain_pos>last_coding_exon[1]) or (acc_gain_pos<first_coding_exon[0]):
			all_flags+=[('OUTSIDE_CODING_REGION',1.0,0.0)]
	else:
		if (acc_gain_pos<last_coding_exon[0]) or (acc_gain_pos>first_coding_exon[1]):
			all_flags+=[('OUTSIDE_CODING_REGION',1.0,0.0)]

	#last exon
	if (last_coding_exon==target_exon) and (first_coding_exon!=last_coding_exon):
		last_coding_frac_lost = (last_coding_portion[1]-last_coding_portion[0]+1)/cds_length
		all_flags+=[('LAST_CODING_EXON',1.0,last_coding_frac_lost)]

	# Cryptic Rescue (with respect to original transcript site)
	nucleotide_difference = np.abs(acc_gain_pos-prior_acceptor_site)

	# inframe cryptic gain, the need for the rescue depends on the strength of the loss prediction
	if ((nucleotide_difference % 3)==0) and (score_table.loc['ACC_LOSS','SCORE']):
		cryptic_rescue_score  = score_table.loc['ACC_LOSS','SCORE']
		if (strand==1):
			if (acc_gain_pos > prior_acceptor_site):
				cryptic_rescue_frac_coding_nuc_impacted= nucleotide_difference/cds_length
			else:
				cryptic_rescue_frac_coding_nuc_impacted=nucleotide_difference/(cds_length+nucleotide_difference)
		else:
			if (acc_gain_pos < prior_acceptor_site):
				cryptic_rescue_frac_coding_nuc_impacted= nucleotide_difference/cds_length
			else:
				cryptic_rescue_frac_coding_nuc_impacted=nucleotide_difference/(cds_length+nucleotide_difference)
		all_flags+=[('CRYPTIC_RESCUE',cryptic_rescue_score,cryptic_rescue_frac_coding_nuc_impacted)]

	return all_flags



def PredictDonorLossRescue(splice_ai_data,variant_pos,strand,exons,utr_5prime,utr_3prime,aa_sequence,met_rescue_threshold=0.9):
	"""
	Given a donor loss, predict rescues as follows:

	1) Original splice site intact: the loss does not impact the normal splice site, or it there is a gain/loss at the same site (normal splice)

	2) Non-coding exon: if the donor impacts a non-coding exon, then it's less likely to have an effect. 

	3) Last exon: a donor upstream of the last exon is more likely to be tolerated. 

	4) Cryptic Rescue: an inframe donor gain rescues a corresponding loss

	5) Exon skipping : An inframe exon is skipped as a result of the loss

	6) Intron retained: an inframe intron can be retained to maintain the reading frame

	7) Methionine resuce: The donor is at the end of the first exon, so that if could be skipped, a methionine in the second could rescue 

	Args:
	    splice_ai_data (dict): Dictionary with the following entries: {'DS_AG','DS_AL','DS_DG','DS_DL','DP_AG','DP_AL','DP_DG':,'DP_DL'}
	    variant_pos (int): genomic position
	    strand (int): 1 or -1
	    exons (list of tuples): transcript exons
	    utr_5prime (list of tuples): tuples defining 5' UTR
	    utr_3prime (list of tuples): tuples defining 3' UTR
	    aa_sequence (str): AA sequence
	    met_rescue_threshold (float, optional): threshold for methionine rescue. If the impacted exon is skipped and is the first exon, opportunity for methionine rescue. At least met_rescue_threshold fraction of the AA sequence must remain. 
	Returns:
	    list of tuples: List of tuples for all possible rescues. Each tuple contains: (label, confidence score, fraction coding nucs impacted)
	"""
	conv_key={'DS_AG':'ACC_GAIN','DS_AL':'ACC_LOSS','DS_DG':'DONOR_GAIN','DS_DL':'DONOR_LOSS','DP_AG':'ACC_GAIN','DP_AL':'ACC_LOSS','DP_DG':'DONOR_GAIN','DP_DL':'DONOR_LOSS'}
	score_table=pd.DataFrame([],columns=['SCORE','GENOME_POS'],index=['ACC_GAIN','ACC_LOSS','DONOR_GAIN','DONOR_LOSS'])
	for key,value in splice_ai_data.items():
		column='SCORE' if key.split('_')[0]=='DS' else 'GENOME_POS'
		if column=='GENOME_POS':
			try:
				genome_pos=variant_pos+value
			except ValueError:
				genome_pos=pd.NA

			score_table.loc[conv_key[key],column]=genome_pos
		else:
			score_table.loc[conv_key[key],column]=value


	donor_loss_pos = score_table.loc['DONOR_LOSS','GENOME_POS']
	donor_gain_pos = score_table.loc['DONOR_GAIN','GENOME_POS']
	downstream_exon=ReturnFirstImpactedExon('DONOR_LOSS',donor_loss_pos,strand,exons)
	target_exon = exons[exons.index(downstream_exon)-1]
	upstream_exon = exons[exons.index(target_exon)-1]
	first_coding_exon,first_coding_portion=ReturnFirstCodingExon(exons,utr_5prime,utr_3prime,strand)
	last_coding_exon,last_coding_portion=ReturnLastCodingExon(exons,utr_5prime,utr_3prime,strand)
	cds_length=ComputeCDSLength(exons,utr_5prime,utr_3prime)

	if (strand==1):
		prior_donor_site = target_exon[1]
	else:
		prior_donor_site=target_exon[0]

	all_flags = []

	# Original splice site intact
	if (prior_donor_site!=donor_loss_pos):
		# if the normal donor site isn't modified, then the only problem is the a new donor site. Since these are delta scores, the intact site score can be approximated by 1 - the gain score
		all_flags += [('ORIG_SPLICESITE_INTACT',1.0-score_table.loc['DONOR_GAIN','SCORE'],0.0)]
	elif score_table.loc['DONOR_GAIN','GENOME_POS']==prior_donor_site:
			all_flags+=[('ORIG_SPLICESITE_INTACT',score_table.loc['DONOR_GAIN','SCORE'],0.0)]

	#outside coding region
	if (strand==1):
		if (prior_donor_site>last_coding_exon[1]) or (prior_donor_site<first_coding_exon[0]):
			all_flags+=[('OUTSIDE_CODING_REGION',1.0,0.0)]
	else:
		if (prior_donor_site<last_coding_exon[0]) or (prior_donor_site>first_coding_exon[1]):
			all_flags+=[('OUTSIDE_CODING_REGION',1.0,0.0)]

	#last exon
	if (last_coding_exon==downstream_exon) and (first_coding_exon!=last_coding_exon):
		last_coding_frac_lost = (last_coding_portion[1]-last_coding_portion[0]+1)/cds_length
		all_flags+=[('LAST_CODING_EXON',1.0,last_coding_frac_lost)]


	# Cryptic Rescue
	nucleotide_difference = np.abs(prior_donor_site-donor_gain_pos)

	if ((nucleotide_difference % 3)==0) and (score_table.loc['DONOR_GAIN','SCORE']>0.0):
		#score for the cryptic rescue 
		cryptic_rescue_score  = score_table.loc['DONOR_GAIN','SCORE']
		if (strand==1):
			if (donor_gain_pos < donor_loss_pos):
				cryptic_rescue_frac_coding_nuc_impacted= nucleotide_difference/cds_length
			else:
				cryptic_rescue_frac_coding_nuc_impacted=nucleotide_difference/(cds_length+nucleotide_difference)
		else:
			if (donor_gain_pos > donor_loss_pos):
				cryptic_rescue_frac_coding_nuc_impacted= nucleotide_difference/cds_length
			else:
				cryptic_rescue_frac_coding_nuc_impacted=nucleotide_difference/(cds_length+nucleotide_difference)
		all_flags+=[('CRYPTIC_RESCUE',cryptic_rescue_score,cryptic_rescue_frac_coding_nuc_impacted)]


	#exon skipping
	if (target_exon!=first_coding_exon) and (target_exon!=last_coding_exon):
		target_exon_size = (target_exon[1]-target_exon[0])+1
	elif (target_exon==first_coding_exon):
		target_exon_size = (first_coding_portion[1]-first_coding_portion[0]+1)
	else:
		target_exon_size = (last_coding_portion[1]-last_coding_portion[0]+1)

	#need to have a meaninful upstream acceptor for this to make sense
	if ((target_exon_size % 3)==0) and (target_exon!=exons[0]):
		if (strand==1):
			upstream_acceptor_position=target_exon[0]
		else:
			upstream_acceptor_position=target_exon[1]
		if score_table.loc['ACC_LOSS','GENOME_POS']==upstream_acceptor_position:
			exon_skip_score = score_table.loc['ACC_LOSS','SCORE']
			exon_skip_frac_coding_nuc_impacted = target_exon_size/cds_length
		else:
			exon_skip_score=0.0
			exon_skip_frac_coding_nuc_impacted=target_exon_size/cds_length
			all_flags+=[('INFRAME_EXON_SKIP',exon_skip_score,exon_skip_frac_coding_nuc_impacted)]


	#intron_retention, need to lose the downstream acceptor. Clearly, this only makes sense if there's a downstream exon. 
	if (target_exon!=exons[-1]):
		if (strand==1):
			intron_size = downstream_exon[0]-target_exon[1]-1
		else:
			intron_size = target_exon[0]-downstream_exon[1]-1
		if (intron_size % 3)==0:
			if (strand==1):
				downstream_acceptor_position=downstream_exon[0]
			else:
				downstream_acceptor_position=downstream_exon[1]

			if score_table.loc['ACC_LOSS','GENOME_POS']==downstream_acceptor_position:
				intron_retention_score = score_table.loc['ACC_LOSS','SCORE']
				intron_retention_frac_coding_nuc_impacted = intron_size/(cds_length+intron_size)
			else:
				intron_retention_score=0.0
				intron_retention_frac_coding_nuc_impacted=intron_size/(cds_length+intron_size)
			all_flags+=[('INFRAME_INTRON_RETENTION',intron_retention_score,intron_retention_frac_coding_nuc_impacted)]

	#finally, possible methionine rescue
	if (target_exon==first_coding_exon):
		first_exon_range = ReturnFirstExonAARange(len(aa_sequence),exons,utr_5prime,utr_3prime,strand)
		if len(aa_sequence[first_exon_range[1]:])/len(aa_sequence) >= met_rescue_threshold:
			#note, off by one indexing is intentional
			new_sequence = aa_sequence[first_exon_range[1]:]
			first_methioninine = new_sequence.find('M')
			new_sequence=new_sequence[new_sequence.find('M'):]
			if len(new_sequence)/len(aa_sequence) >= met_rescue_threshold:
				possible_met_rescue=True
				met_rescue_nuc_impacted = 1.0-(len(aa_sequence[first_exon_range[1]:])/len(aa_sequence))
				all_flags+=[('POSSIBLE_MET_RESCUE',1.0,met_rescue_nuc_impacted)]
	return all_flags




def PredictDonorGainRescue(splice_ai_data,variant_pos,strand,exons,utr_5prime,utr_3prime,met_rescue_threshold=0.9):
	"""
	Given a donor gain event, assess the following rescue possibilities:

	1) Original splice site intact: the gain simply rescues/impacts the same normal splice site
	2) Non-coding exon: if the donor impacts a non-coding region of the transcript, then it's unlikely to have an effect. 
	3) Last exon: a donor gain upstream of the last exon is more likely to be tolerated. 
	4) Cryptic Rescue: the gain rescues an inframe loss

 	Becuase the gain of a donor is inducing a new splice site, its less likely that exon skipping, intron retention, or methionine restart could rescue. 
	
	Args:
	    splice_ai_data (dict): Dictionary with the following entries: {'DS_AG','DS_AL','DS_DG','DS_DL','DP_AG','DP_AL','DP_DG':,'DP_DL'}
	    variant_pos (int): genomic position
	    strand (int): 1 or -1
	    exons (list of tuples): transcript exons
	    utr_5prime (list of tuples): tuples defining 5' UTR
	    utr_3prime (list of tuples): tuples defining 3' UTR
	
	Returns:
	    list of tuples: List of tuples for all possible rescues. Each tuple contains: (label, confidence score, fraction coding nucs impacted)
	"""
	conv_key={'DS_AG':'ACC_GAIN','DS_AL':'ACC_LOSS','DS_DG':'DONOR_GAIN','DS_DL':'DONOR_LOSS','DP_AG':'ACC_GAIN','DP_AL':'ACC_LOSS','DP_DG':'DONOR_GAIN','DP_DL':'DONOR_LOSS'}
	score_table=pd.DataFrame([],columns=['SCORE','GENOME_POS'],index=['ACC_GAIN','ACC_LOSS','DONOR_GAIN','DONOR_LOSS'])
	for key,value in splice_ai_data.items():
		column='SCORE' if key.split('_')[0]=='DS' else 'GENOME_POS'
		if column=='GENOME_POS':
			try:
				genome_pos=variant_pos+value
			except ValueError:
				genome_pos=pd.NA

			score_table.loc[conv_key[key],column]=genome_pos
		else:
			score_table.loc[conv_key[key],column]=value

	donor_loss_pos = score_table.loc['DONOR_LOSS','GENOME_POS']
	donor_gain_pos = score_table.loc['DONOR_GAIN','GENOME_POS']
	downstream_exon=ReturnFirstImpactedExon('DONOR_GAIN',donor_gain_pos,strand,exons)
	target_exon = exons[exons.index(downstream_exon)-1]
	upstream_exon = exons[exons.index(target_exon)-1]
	first_coding_exon,first_coding_portion=ReturnFirstCodingExon(exons,utr_5prime,utr_3prime,strand)
	last_coding_exon,last_coding_portion=ReturnLastCodingExon(exons,utr_5prime,utr_3prime,strand)
	cds_length=ComputeCDSLength(exons,utr_5prime,utr_3prime)


	if (strand==1):
		prior_donor_site = target_exon[1]
	else:
		prior_donor_site=target_exon[0]

	all_flags = []

	# Original splice site intact
	if prior_donor_site==donor_gain_pos:
		# the gain is at the target site. Just return the donor gain score (mirrors donor loss issue)
		 all_flags+=[('ORIG_SPLICESITE_INTACT',score_table.loc['DONOR_GAIN','SCORE'],0.0)]

	elif (prior_donor_site!=score_table.loc['DONOR_LOSS','GENOME_POS']):
		# if the normal donor site isn't modified, then the only problem is the a new donor site. Since these are delta scores, the intact site score can be approximated by 1 - the gain score
		 all_flags+=[('ORIG_SPLICESITE_INTACT',1.0-score_table.loc['DONOR_GAIN','SCORE'],0.0)]

	#outside coding region
	if (strand==1):
		if (prior_donor_site>last_coding_exon[1]) or (prior_donor_site<first_coding_exon[0]):
			all_flags+=[('OUTSIDE_CODING_REGION',1.0,0.0)]
	else:
		if (prior_donor_site<last_coding_exon[0]) or (prior_donor_site>first_coding_exon[1]):
			all_flags+=[('OUTSIDE_CODING_REGION',1.0,0.0)]

	#last exon
	if (last_coding_exon==downstream_exon) and (first_coding_exon!=last_coding_exon):
		last_coding_frac_lost = (last_coding_portion[1]-last_coding_portion[0]+1)/cds_length
		all_flags+=[('LAST_CODING_EXON',1.0,last_coding_frac_lost)]


	nucleotide_difference = np.abs(donor_loss_pos-prior_donor_site)

	# inframe cryptic gain
	if ((nucleotide_difference % 3)==0) and (score_table.loc['DONOR_LOSS','SCORE']>0.0):
		#the cryptic rescue role for a splice gain is proportional to the score of the loss at the old site.
		cryptic_rescue_score  = score_table.loc['DONOR_LOSS','SCORE']

		if (strand==1):
			if (donor_gain_pos < donor_loss_pos):
				cryptic_rescue_frac_coding_nuc_impacted= nucleotide_difference/cds_length
			else:
				cryptic_rescue_frac_coding_nuc_impacted=nucleotide_difference/(cds_length+nucleotide_difference)
		else:
			if (donor_gain_pos > donor_loss_pos):
				cryptic_rescue_frac_coding_nuc_impacted= nucleotide_difference/cds_length
			else:
				cryptic_rescue_frac_coding_nuc_impacted=nucleotide_difference/(cds_length+nucleotide_difference)
		all_flags+=[('CRYPTIC_RESCUE',cryptic_rescue_score,cryptic_rescue_frac_coding_nuc_impacted)]

	return all_flags



def PredictEscapeNMD(variant_pos,exon_list,utr_5prime,utr_3prime,strand):
	"""Based on the decision tree provided in PMID: 31659324, Figure 1c. 
	
	Args:
	    variant_pos (int): genomic positon of the variant

	    exons (list of tuples): transcript exons
	    utr_5prime (list of tuples): tuples defining 5' UTR
	    utr_3prime (list of tuples): tuples defining 3' UTR
	    strand (int): 1 or -1
	Returns:
	    (bool,str): Returns the whether the transcript is predicted to escape NMD, along with a flag indicating the reason. 
	"""
	last_coding_exon,last_coding_exon_trimmed=ReturnLastCodingExon(exon_list,utr_5prime,utr_3prime,strand)
	first_coding_exon,first_coding_exon_trimmed=ReturnFirstCodingExon(exon_list,utr_5prime,utr_3prime,strand)
	if ((variant_pos>=last_coding_exon_trimmed[0]) and (variant_pos<=last_coding_exon_trimmed[1])) and (last_coding_exon!=first_coding_exon):
		#last exon
		nmd_escape=True
		flag='LAST_CODING_EXON'

	else:
		if strand==1:
			distance_from_coding_start=(variant_pos-first_coding_exon_trimmed[0])
		else:
			distance_from_coding_start=(first_coding_exon_trimmed[1]-variant_pos)

		assert distance_from_coding_start>=0,"Variant {0:s} has negative distance from Start".format(variant_pos)

		#in the first exon within 150bp of coding start
		if ((variant_pos>=first_coding_exon_trimmed[0]) and (variant_pos<=first_coding_exon_trimmed[1])) and (distance_from_coding_start<=150):
			nmd_escape=True
			flag='FIRST_EXON_LEQ_150NT_FROM_START'
		else:
			#compute exon length
			for containing_exon in exon_list:
				if (variant_pos>=containing_exon[0]) and (variant_pos<=containing_exon[1]):
					break
			containing_exon_number=exon_list.index(containing_exon)+1
			containing_exon_length = (containing_exon[1]-containing_exon[0])+1
			if containing_exon_length>407:
				nmd_escape=True
				flag='LARGE_EXON'
			else:
				penultimate_coding_exon = exon_list.index(last_coding_exon)

				if strand==1:
					last_exon_junction=exon_list[penultimate_coding_exon-1][1]
					dist_from_last_ej = (last_exon_junction-variant_pos)
				else:
					last_exon_junction=exon_list[penultimate_coding_exon-1][0]
					dist_from_last_ej = (variant_pos-last_exon_junction)
				assert dist_from_last_ej>=0,"Variant {0:s} has negative distance from last EJ".format(variant_id)
				if (containing_exon_number==penultimate_coding_exon) and (dist_from_last_ej<=50):
					nmd_escape=True
					flag='LEQ_50NT_FROM_LAST_EJ'
				else:
					nmd_escape=False
					flag=pd.NA

	return nmd_escape,flag



def PredictAAImpacted_NoSplice(consequence,position_info,ref_residues,alt_residues,aa_sequence,exon_list,utr_5prime,utr_3prime,strand,met_rescue_threshold=0.9):
	"""Predicts the number/fraction of AA impacted assuming that the transcript escapes NMD. Used only for FRAMESHIFT and stop gain
	
	Args:
	    consequence (str): Description
	    position_info (str): AA position info. Can be a single value or a range.
	    ref_residues (str): ref aa residue(s)
	    alt_residues (str): alt aa residue(s)
	    aa_sequence (str): full AA sequence
	    exons (list of tuples): transcript exons
	    utr_5prime (list of tuples): tuples defining 5' UTR
	    utr_3prime (list of tuples): tuples defining 3' UTR
	    sstrand (int): 1 or -1
	    met_rescue_threshold (float, optional): threshold for methionine rescue. If the impacted exon is the first exon, opportunity for methionine rescue (restricted to this exon only). At least met_rescue_threshold fraction of the AA sequence must remain. 
	
	Returns:
	    tuple: num_aa_lost,frac_aa_lost,flag
	

	"""
	last_exon_AA_range=ReturnLastExonAARange(len(aa_sequence),exon_list,utr_5prime,utr_3prime,strand)
	first_exon_AA_range=ReturnFirstExonAARange(len(aa_sequence),exon_list,utr_5prime,utr_3prime,strand)

	#identify which residues are impacted
	if consequence=='STOP_GAINED':
		if len(position_info.split('-'))==1:
			#single AA, easy
			residue_pos=int(position_info)
			residue=ref_residues

		else:
			#multi-residue AA 
			res_start=int(position_info.split('-')[0])
			res_stop=int(position_info.split('-')[1])
			pos_sequence=list(range(res_start,res_stop+1))
			reference_residues = ref_residues.replace('-','')
			mutated_residues = alt_residues
			try:
				residue_pos=pos_sequence[mutated_residues.index('*')]
				residue=reference_residues[mutated_residues.index('*')]
			except IndexError:
				#is indel, default to lowest normal residue
				assert reference_residues=='',"Indel error for variant {0:s}.".format(variant_id)
				residue_pos=res_start+mutated_residues.index('*')
				residue = aa_sequence[residue_pos-1]

		# residue position lies outside the translated AA sequence. This can happen when there are discrepancies between AA sequence listed in the databases vs what is actually predicted based on the cDNA sequence. In such cases, defaults to TX_TRANSLATION_ERROR, which are likely to be tolerated. 
		if residue_pos>len(aa_sequence):
			num_aa_lost=0
			frac_aa_lost=0.0
			flag='TX_TRANSLATION_ERROR'
		else:
			#again, if there's a discrepancy between the AA provided by VEP and the one predicted from the sequence, there is likely a TX_TRANSLATION_ERROR. These are rare but do occur. 
			if aa_sequence[residue_pos-1]!=residue:
				num_aa_lost=0
				frac_aa_lost=0.0
				flag='TX_TRANSLATION_ERROR'
			else:
				#everything is lost downstream of the stop codon. Note the off-by-1 indexing to account AA positions starting at 1
				num_aa_lost=len(aa_sequence[residue_pos-1:])
				frac_aa_lost=num_aa_lost/len(aa_sequence)


				#now the flags
				#if the mutation occurs in the first exon, has a downstream MET in this exon, and the total loss is less than met_rescue_threshold% of the coding sequence. Again, watch the off-by-1 indexing carefully. 
				if (residue_pos<=first_exon_AA_range[1]) and ('M' in str(aa_sequence[residue_pos:first_exon_AA_range[1]])):
					new_met_rel_position = str(aa_sequence[residue_pos:first_exon_AA_range[1]]).index('M')
					new_met_abs_position = residue_pos+new_met_rel_position+1
					#just a check to make sure everything is internally consistent. 
					assert aa_sequence[new_met_abs_position-1]=='M',"Possible new methionine start does not match M in sequence."
					if len(aa_sequence[new_met_abs_position-1:])/len(aa_sequence)>=met_rescue_threshold:
						flag='POSSIBLE_MET_RESCUE'
					else:
						flag=pd.NA

				#if in the last exon and there are more than 1 exons; note this is redundant with NMD escape prediction, but include for symmetry with frameshifts
				elif (last_exon_AA_range[0]<=residue_pos) and (last_exon_AA_range[1]>=residue_pos) and (last_exon_AA_range[0]>first_exon_AA_range[1]):
					flag='LAST_CODING_EXON'
				else:
					flag=pd.NA

	elif consequence=='FRAMESHIFT':
		aa_start=int(position_info.split('-')[0])
		altered_aa=alt_residues
		assert (('X' in altered_aa) or ('*' in altered_aa)) or (altered_aa==''),"Frameshift variant with ID {0:s} does not contain expected altered AA sequence".format(variant_id)
		pos_sequence=list(range(aa_start,aa_start+len(altered_aa)+1))

		#find induced stop codon by the frameshift. If not present, then default to the first AA in the altered sequence. 
		if '*' in altered_aa:
			affected_pos=pos_sequence[altered_aa.index('*')]	
		elif 'X' in altered_aa:
			affected_pos=pos_sequence[altered_aa.index('X')]
		else:
			affected_pos=aa_start-1 #use the last normal AA as defualt if no additional information

		num_aa_lost=len(aa_sequence[affected_pos-1:])
		frac_aa_lost=num_aa_lost/len(aa_sequence)

		#methionine rescue, just like for stop codons, is theoretically possible 
		if (affected_pos<=first_exon_AA_range[1]) and ('M' in str(aa_sequence[affected_pos:first_exon_AA_range[1]])):
			new_met_rel_position = str(aa_sequence[affected_pos:first_exon_AA_range[1]]).index('M')
			new_met_abs_position = affected_pos+new_met_rel_position+1
			assert aa_sequence[new_met_abs_position-1]=='M',"Possible new methionine start does not match M in sequence."
			if len(aa_sequence[new_met_abs_position-1:])/len(aa_sequence)>=met_rescue_threshold:
				flag='POSSIBLE_MET_RESCUE'
			else:
				flag=pd.NA
		#frameshifts in the last exon are more likely to be tolerted, either due to NMD escape or just because only the last few AAs are altered. This could of course be made a lot more precise. 
		elif (last_exon_AA_range[0]<=affected_pos) and (last_exon_AA_range[1]>=affected_pos) and (last_exon_AA_range[0]>first_exon_AA_range[1]):
			flag='LAST_CODING_EXON'
		else:
			flag=pd.NA
	else:
		raise ValueError('Consequence {0:s} not recognized'.format(consequence))
	return num_aa_lost,frac_aa_lost,flag

