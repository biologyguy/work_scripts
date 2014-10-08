#!/usr/bin/python

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import random
import math


NUM_IN_GEN = 10000
NUM_GENS = 100    
ID_PADDING = 7
GEN_PADDING = 3


def test_mutation(residue):
    
    blosum62_matrix = {'A': {'A': '0.917974003', 'C': '0.932317346', 'E': '0.939489018', 'D': '0.925145675', 'G': '0.94666069', 'F': '0.96458987', 'I': '0.952039444', 'H': '0.948453608', 'K': '0.959211116', '*': '1', 'M': '0.962796952', 'L': '0.95562528', 'N': '0.923352757', 'Q': '0.935903182', 'P': '0.968175706', 'S': '0.98251905', 'R': '0.921559839', 'T': '0.989690722', 'W': '0.990587181', 'V': '0.999551771', 'Y': '0.992380099'}, 'C': {'A': '0.000243791', 'C': '0.998902941', 'E': '0.998948652', 'D': '0.000335213', 'G': '0.998979125', 'F': '0.999466707', 'I': '0.999131495', 'H': '0.999009599', 'K': '0.999283864', '*': '1', 'M': '0.99940576', 'L': '0.99925339', 'N': '0.000304739', 'Q': '0.998933415', 'P': '0.999497181', 'S': '0.999619077', 'R': '0.000274265', 'T': '0.999740972', 'W': '0.99980192', 'V': '0.999984763', 'Y': '0.999862868'}, 'E': {'A': '0.001830664', 'C': '0.02402746', 'E': '0.97597254', 'D': '0.023798627', 'G': '0.976887872', 'F': '0.990160183', 'I': '0.981006865', 'H': '0.980549199', 'K': '0.988787185', '*': '1', 'M': '0.989702517', 'L': '0.981464531', 'N': '0.009153318', 'Q': '0.038672769', 'P': '0.991990847', 'S': '0.995652174', 'R': '0.005491991', 'T': '0.997482838', 'W': '0.997940503', 'V': '0.999771167', 'Y': '0.998855835'}, 'D': {'A': '0.000477156', 'C': '0.982225933', 'E': '0.991769056', 'D': '0.981987355', 'G': '0.992723369', 'F': '0.995467017', 'I': '0.993916259', 'H': '0.993677681', 'K': '0.99498986', '*': '1', 'M': '0.995228439', 'L': '0.994035548', 'N': '0.004771561', 'Q': '0.984134558', 'P': '0.996421329', 'S': '0.998329953', 'R': '0.000954312', 'T': '0.999284266', 'W': '0.999403555', 'V': '0.999880711', 'Y': '0.999642133'}, 'G': {'A': '0.001929338', 'C': '0.005546847', 'E': '0.006511516', 'D': '0.005305679', 'G': '0.99433257', 'F': '0.99602074', 'I': '0.994935488', 'H': '0.994814904', 'K': '0.995538406', '*': '1', 'M': '0.995779573', 'L': '0.995056071', 'N': '0.00434101', 'Q': '0.006029181', 'P': '0.996503075', 'S': '0.998432413', 'R': '0.002411672', 'T': '0.998914747', 'W': '0.999397082', 'V': '0.999879416', 'Y': '0.999638249'}, 'F': {'A': '0.000473485', 'C': '0.001657197', 'E': '0.002130682', 'D': '0.001183712', 'G': '0.002367424', 'F': '0.978929924', 'I': '0.005208333', 'H': '0.003314394', 'K': '0.007339015', '*': '1', 'M': '0.009232955', 'L': '0.007102273', 'N': '0.00094697', 'Q': '0.001893939', 'P': '0.979048295', 'S': '0.97952178', 'R': '0.000710227', 'T': '0.979995265', 'W': '0.983783144', 'V': '0.999881629', 'Y': '0.998934659'}, 'I': {'A': '0.003412969', 'C': '0.009385666', 'E': '0.01109215', 'D': '0.005972696', 'G': '0.011518771', 'F': '0.934726962', 'I': '0.88609215', 'H': '0.012372014', 'K': '0.914249147', '*': '1', 'M': '0.927901024', 'L': '0.913395904', 'N': '0.005119454', 'Q': '0.010238908', 'P': '0.935580205', 'S': '0.937286689', 'R': '0.004266212', 'T': '0.940699659', 'W': '0.941552901', 'V': '0.999573379', 'Y': '0.94496587'}, 'H': {'A': '0.000121297', 'C': '0.001880098', 'E': '0.002850472', 'D': '0.00181945', 'G': '0.002971768', 'F': '0.997361798', 'I': '0.996694666', 'H': '0.996634018', 'K': '0.996997908', '*': '1', 'M': '0.997119204', 'L': '0.996755314', 'N': '0.001576857', 'Q': '0.002365285', 'P': '0.997483094', 'S': '0.997725688', 'R': '0.000606483', 'T': '0.997846984', 'W': '0.997968281', 'V': '0.999969676', 'Y': '0.999909028'}, 'K': {'A': '0.001848002', 'C': '0.022638023', 'E': '0.037422037', 'D': '0.022176022', 'G': '0.038346038', 'F': '0.99006699', 'I': '0.040656041', 'H': '0.04019404', 'K': '0.987756988', '*': '1', 'M': '0.98960499', 'L': '0.041580042', 'N': '0.02032802', 'Q': '0.03003003', 'P': '0.991914992', 'S': '0.995610996', 'R': '0.016632017', 'T': '0.997458997', 'W': '0.997920998', 'V': '0.999769', 'Y': '0.998844999'}, '*': {}, 'M': {'A': '0.001842893', 'C': '0.00691085', 'E': '0.011518083', 'D': '0.005067957', 'G': '0.011978807', 'F': '0.984105045', 'I': '0.020271827', 'H': '0.012900253', 'K': '0.036857867', '*': '1', 'M': '0.980419258', 'L': '0.035014974', 'N': '0.004607233', 'Q': '0.010596637', 'P': '0.985026492', 'S': '0.986869385', 'R': '0.003685787', 'T': '0.988712278', 'W': '0.990555172', 'V': '0.999769638', 'Y': '0.992398065'}, 'L': {'A': '0.003498032', 'C': '0.010056843', 'E': '0.012680367', 'D': '0.006558811', 'G': '0.013117621', 'F': '0.974202011', 'I': '0.041976388', 'H': '0.013992129', 'K': '0.939221688', '*': '1', 'M': '0.967205947', 'L': '0.937472672', 'N': '0.006121557', 'Q': '0.011805859', 'P': '0.975076519', 'S': '0.976825536', 'R': '0.005247049', 'T': '0.980323568', 'W': '0.982072584', 'V': '0.999562746', 'Y': '0.985570617'}, 'N': {'A': '0.000475511', 'C': '0.980266286', 'E': '0.984070376', 'D': '0.980028531', 'G': '0.98597242', 'F': '0.992867332', 'I': '0.990014265', 'H': '0.98977651', 'K': '0.992154066', '*': '1', 'M': '0.992629577', 'L': '0.990252021', 'N': '0.976224441', 'Q': '0.982168331', 'P': '0.993342844', 'S': '0.997146933', 'R': '0.002377556', 'T': '0.999048978', 'W': '0.999167855', 'V': '0.999881122', 'Y': '0.999643367'}, 'Q': {'A': '0.001835283', 'C': '0.016976371', 'E': '0.971323698', 'D': '0.01651755', 'G': '0.97224134', 'F': '0.98875889', 'I': '0.976370727', 'H': '0.975911906', 'K': '0.984629502', '*': '1', 'M': '0.988300069', 'L': '0.977288369', 'N': '0.012846983', 'Q': '0.956641432', 'P': '0.990594173', 'S': '0.99426474', 'R': '0.009176417', 'T': '0.996100023', 'W': '0.997017665', 'V': '0.99977059', 'Y': '0.998852948'}, 'P': {'A': '0.000485584', 'C': '0.001578149', 'E': '0.002549317', 'D': '0.001456753', 'G': '0.002792109', 'F': '0.004066768', 'I': '0.003156297', 'H': '0.003034901', 'K': '0.003763278', '*': '1', 'M': '0.00400607', 'L': '0.003277693', 'N': '0.000971168', 'Q': '0.002063733', 'P': '0.998543247', 'S': '0.999028832', 'R': '0.000728376', 'T': '0.999514416', 'W': '0.999575114', 'V': '0.999939302', 'Y': '0.99969651'}, 'S': {'A': '0.013992129', 'C': '0.041976388', 'E': '0.055968518', 'D': '0.038478356', 'G': '0.062964582', 'F': '0.08220376', 'I': '0.068211631', 'H': '0.066462615', 'K': '0.076956712', '*': '1', 'M': '0.080454744', 'L': '0.069960647', 'N': '0.031482291', 'Q': '0.048972453', 'P': '0.085701793', 'S': '0.981198076', 'R': '0.017490162', 'T': '0.995190206', 'W': '0.996064714', 'V': '0.999562746', 'Y': '0.99781373'}, 'R': {'A': '0.001859168', 'C': '0.958865908', 'E': '0.970020916', 'D': '0.958401116', 'G': '0.9709505', 'F': '0.993260516', 'I': '0.975133628', 'H': '0.974668836', 'K': '0.990936556', '*': '1', 'M': '0.992795724', 'L': '0.976063212', 'N': '0.957471531', 'Q': '0.96630258', 'P': '0.9941901', 'S': '0.996049268', 'R': '0.953753195', 'T': '0.997908436', 'W': '0.998373228', 'V': '0.999767604', 'Y': '0.999302812'}, 'T': {'A': '0.00374094', 'C': '0.01309329', 'E': '0.01683423', 'D': '0.01122282', 'G': '0.017769465', 'F': '0.027121814', 'I': '0.02057517', 'H': '0.0187047', 'K': '0.024316109', '*': '1', 'M': '0.026186579', 'L': '0.022445639', 'N': '0.00935235', 'Q': '0.01496376', 'P': '0.028992284', 'S': '0.036474164', 'R': '0.00561141', 'T': '0.994154781', 'W': '0.995090016', 'V': '0.999766191', 'Y': '0.996025251'}, 'W': {'A': '7.62515E-06', 'C': '3.81257E-05', 'E': '6.10012E-05', 'D': '2.28754E-05', 'G': '7.62515E-05', 'F': '0.000274505', 'I': '9.91269E-05', 'H': '9.15018E-05', 'K': '0.000122002', '*': '1', 'M': '0.000152503', 'L': '0.000114377', 'N': '1.90629E-05', 'Q': '5.3376E-05', 'P': '0.000278318', 'S': '0.000285943', 'R': '1.52503E-05', 'T': '0.000301193', 'W': '0.999744558', 'V': '0.999996187', 'Y': '0.999988562'}, 'V': {'A': '0.006869901', 'C': '0.012881065', 'E': '0.016316015', 'D': '0.009446114', 'G': '0.017174753', 'F': '0.105624732', 'I': '0.072992701', 'H': '0.018033491', 'K': '0.088449979', '*': '1', 'M': '0.102189781', 'L': '0.086732503', 'N': '0.008587377', 'Q': '0.01459854', 'P': '0.107342207', 'S': '0.109059682', 'R': '0.007728639', 'T': '0.115929584', 'W': '0.116788321', 'V': '0.999570631', 'Y': '0.120223272'}, 'Y': {'A': '0.000239249', 'C': '0.001076619', 'E': '0.001794366', 'D': '0.000837371', 'G': '0.00191399', 'F': '0.015072672', 'I': '0.006220468', 'H': '0.00574197', 'K': '0.006938214', '*': '1', 'M': '0.007416712', 'L': '0.006698965', 'N': '0.000717746', 'Q': '0.001555117', 'P': '0.015192296', 'S': '0.015431545', 'R': '0.000478498', 'T': '0.015670794', 'W': '0.019498774', 'V': '0.999940188', 'Y': '0.99946169'}}

    order = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "*"]
    
    rand_num = random.random()
    for aa in order:
        temp = blosum62_matrix[residue]
        if rand_num < float(temp[aa]):
            if aa == "*":
                indel_size = int(math.floor(1 / random.random()))
                return indel_size
            
            return aa
    

def calc_fitness(seq):  # (SeqRecord object)
    #Parameters: 4 transmembrane domains, C residues in ECLs, (Y|W)YQW in 2nd TM domain, FP(K|R)(M|V|L|I) in 2nd ECL, Methionine at start
    fitness = 100
    
    #big fitness loss if Methionine is not first residue
    if seq.seq[0] != "M":
        fitness *= 0.25
    
    #Select against large size descrepencies
    seq_length = len(seq.seq)
    if seq_length > 550 or seq_length < 350:
        fitness *= 0.5
    
    return fitness


def mutate(seq, pos_list, aa_position_count):
    new_seq = ""
    amino_acids = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
    deletion_count = 0
    pos_list_counter = -1
    positions = ""
    
    for aa in seq:
        pos_list_counter += 1
        if deletion_count > 0:
            deletion_count -= 1
            continue
            
        aa = test_mutation(aa)
        
        if type(aa) == str:
            new_seq += aa
            positions += pos_list[pos_list_counter] + ","
            
        else:
            if random.random() < 0.5:
                #Deletion
                deletion_count = aa
            
            else:
                #Insertion
                insertion = ""
                for x in range(0, aa):
                    insertion += random.choice(amino_acids)
                    positions += str(aa_position_count) + ","
                    aa_position_count += 1
                new_seq += insertion

    output_dic = {"seq": new_seq, "pos": positions, "pos_count": aa_position_count}
    return output_dic
    

sequence_id = 0
gen_id = 0

#Generation 0
ancestor = Seq("MLLLGSLGTIKNLSIFKDLSLDDWLDQMNRTFMFLLLCFMGTIVAVSQYTGKNISCDGFTKFGEDFSQDYCWTQGLYTIKEAYDLPESQIPYPGIIPENVPACREHALKNGGKIVCPPEDQVKPLTRARHLWYQWIPFYFWVIAPVFYLPYMFVKRMGLDRMKPLLKIMSDYYHCTTETPSEEIIVKCADWVYNSIVDRLSEGSSWTSWRNRHGLGLAVLVSKFMYLGGSVLVMMMTTLMFQVGDFKTYGIEWLRQFPNPENYSTSVKHKLFPKMVACEIKRWGTTGLEEENGMCVLAPNVIYQYIFLIMWFALAITICTNFGNIFFYLFKLTATRYTYNKLVATGHFSHKHPGWKFMYYRIGTSGRVLLNIVAQNTNPIIFGAIMEKLTPSVIKHLRIGHVPGEYLTDPA", 'generic_protein')

fasta_output = open("seq_files/gen000.fasta", 'w')
aa_positions_output = open("seq_files/gen000.pos", 'w')
aa_position_count = len(ancestor)


for _ in range(0, NUM_IN_GEN):
    next_id = str(sequence_id)
    sequence_id += 1
    sequence = SeqRecord(ancestor, id=next_id.zfill(ID_PADDING), description="parent 0000")
    SeqIO.write(sequence, fasta_output, "fasta")
    aa_positions_output.write(next_id.zfill(ID_PADDING) + "\n")
    for i in range(0, aa_position_count):
        aa_positions_output.write(str(i) + ",")
    aa_positions_output.write("\n")

fasta_output.close()
aa_positions_output.close()

#Offspring generations
for _ in range(0, NUM_GENS):
    #need to open the parental positions file before iterating gen_id
    next_gen = str(gen_id)
    parental_positions_file = open("seq_files/gen" + next_gen.zfill(GEN_PADDING) + ".pos", 'r')
    gen_id += 1
    
    parental_input = SeqIO.parse("seq_files/gen" + next_gen.zfill(GEN_PADDING) + ".fasta", "fasta")
    parental_input_dict = {}
    sum_fitness = 0.0
    fitness_dict = {}
    
    #get fitness values for each sequence
    for next_seq in parental_input:
        parental_input_dict[next_seq.id] = next_seq
        fitness = calc_fitness(next_seq)
        fitness_dict[next_seq.id] = fitness
        sum_fitness += fitness
    
    #Change fitness values from ints to cumulative probabilities
    cumu_sum_prob = 0.0
    for next_seq_id in fitness_dict:
        fitness_prob = fitness_dict[next_seq_id] / sum_fitness
        fitness_dict[next_seq_id] = fitness_prob + cumu_sum_prob
        cumu_sum_prob += fitness_prob
    
    #Build next generation
    next_gen = str(gen_id)    
    fasta_output = open("seq_files/gen" + next_gen.zfill(GEN_PADDING) + ".fasta", "w")
    aa_positions_output = open("seq_files/gen" + next_gen.zfill(GEN_PADDING) + ".pos", 'w')
    
    parental_positions_dict = {}
    for next_line in parental_positions_file:
        if len(next_line.strip()) == ID_PADDING:
            parental_positions_dict[next_line.strip()] = []
            dict_id = next_line.strip()
        else:
            positions = next_line.strip().split(',')
            positions.pop(-1)
             
            for next_pos in positions:
                parental_positions_dict[dict_id].append(next_pos)

    #Start by selecting unmutated parental sequences
    for _ in range(0, int(NUM_IN_GEN * 0.1)):
        rand_num = random.random()
        
        for next_seq_id in fitness_dict:
            if rand_num < float(fitness_dict[next_seq_id]):
                next_id = str(sequence_id)
                sequence_id += 1
                child = SeqRecord(parental_input_dict[next_seq_id].seq, id=next_id.zfill(ID_PADDING), description="parent " + parental_input_dict[next_seq_id].id)
                SeqIO.write(child, fasta_output, "fasta")

                aa_positions_output.write(child.id + "\n")
                for pos_key in parental_positions_dict[next_seq_id]:
                    aa_positions_output.write(pos_key + ",")
                
                aa_positions_output.write("\n")
                break
        
    #Now add remaining sequences with mutations
    for _ in range(0, int(NUM_IN_GEN * 0.9)):
        rand_num = random.random()
        
        for next_seq_id in fitness_dict:
            if rand_num < float(fitness_dict[next_seq_id]):
                next_id = str(sequence_id)
                sequence_id += 1
                
                mutation = mutate(parental_input_dict[next_seq_id].seq, parental_positions_dict[parental_input_dict[next_seq_id].id], aa_position_count)
                aa_position_count = mutation['pos_count']
                mutated_seq = Seq(mutation["seq"], "generic_protein")
                
                child = SeqRecord(mutated_seq, id=next_id.zfill(ID_PADDING), description="parent " + parental_input_dict[next_seq_id].id)
                SeqIO.write(child, fasta_output, "fasta")
                
                aa_positions_output.write(child.id + "\n" + mutation['pos'] + "\n")
                
                break

    fasta_output.close()
    aa_positions_output.close()
    parental_positions_file.close()

print("Done")


'''
#Use this to get the blosum matrix
blosum = open("BLOSUM62_better","r")
order = ["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V","*"]
final = {"A":{},"R":{},"N":{},"D":{},"C":{},"Q":{},"E":{},"G":{},"H":{},"I":{},"L":{},"K":{},"M":{},"F":{},"P":{},"S":{},"T":{},"W":{},"Y":{},"V":{},"*":{}}

for next_line in blosum:
    clean_row = next_line.strip()
    values = clean_row.split("\t")
    print values
    current_aa = values.pop(0)
    
    counter = 0    
    for value in values:
        final[current_aa][order[counter]] = value
        counter += 1
    
print final
'''