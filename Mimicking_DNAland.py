# -*- coding: utf-8 -*-

def relu(x):
    s = np.where(x < 0, 0, x)
    return s

def psigmoid(x):
    s = 1-4**(-x**2)
    return s
def transpose(matrix):
    return zip(*matrix)
import numpy as np
import gzip
genome_dict = {}
with open(r"H:\BAM_Analysis_Kit_1_8\BAM_Analysis_Kit_1_8\DA246\filtered-autosomal-o37-results.csv\dnl112247_tpy.txt", 'r') as idx_f:
    for line in idx_f:
        if (line.startswith('#') or line.startswith('\n')) is not True:
            fields = line.strip().split('\t')
            # index_pos = fields[0]
            rsid = fields[0]
            chromosome = fields[1]
            position = fields[2]
            genotype = fields[3]
            # start_pos = int(index_pos) * 2
            ch_rs_name = chromosome + ':' + position
            genome_dict[ch_rs_name] = {
                'chromosome': chromosome,
                'position': position,
                'genotype': genotype,
            }
user_genome = genome_dict
dnaland_tem = gzip.open(r'refpanel_public_ref.gz')
nation_number_total = [0] * 50
one_line_total_number_log = [0] * 50
unweighted_nation_number = [0] * 50
matchNO = 0
for dnaland_tem_line in dnaland_tem:
    dnaland_tem_line_fields = dnaland_tem_line.strip().split( )#the space
    dnaland_tem_line_fields_ch_rs_name = dnaland_tem_line_fields[1]
    dnaland_tem_line_fields_A1 = dnaland_tem_line_fields[4]
    dnaland_tem_line_fields_A2 = dnaland_tem_line_fields[5]
    if dnaland_tem_line_fields_ch_rs_name in user_genome:
        user_genotype = user_genome[dnaland_tem_line_fields_ch_rs_name]['genotype']
        one_line_total_number = [0] * 42
        nation_frequence = [0] * 42
        for nation_number in range(6, 48):
            nation_number_tem = 0
            dnaland_tem_line_fields_nation = dnaland_tem_line_fields[nation_number]
            dnaland_tem_line_fields_nation_number_fields = dnaland_tem_line_fields_nation.strip().split(',')
            nation_number_A1 = dnaland_tem_line_fields_nation_number_fields[0]
            nation_number_A2 = dnaland_tem_line_fields_nation_number_fields[1]
            A1_factor = 0.000000001
            A2_factor = 0.000000001
            for dnabase in user_genotype:
                if dnabase == dnaland_tem_line_fields_A1:
                    A1_factor = A1_factor + 1
                if dnabase == dnaland_tem_line_fields_A2:
                    A2_factor = A2_factor + 1
                if not (nation_number_A1 + nation_number_A2) < 0.1:
                    nation_number_like = ((float(nation_number_A1)+0.00000001)*A1_factor) + ((float(nation_number_A2)+0.00000001) *A2_factor)
                    one_line_total_number[nation_number-6] = nation_number_like
                    nation_frequence[nation_number - 6] = float(nation_number_A1) +float(nation_number_A2)
        nation_number_total = map(lambda (a, b): a + b, zip(nation_number_total, one_line_total_number))
        unweighted_nation_number = map(lambda (a, b): a + b, zip(unweighted_nation_number, nation_frequence))

total_unweighted_nation_number = sum(unweighted_nation_number)
#nation_number_total_per = [c/total_number for c in nation_number_total]
per_unweighted_nation_number = np.multiply(nation_number_total, (1/total_unweighted_nation_number))
nation_name = ['东北亚NEASIA','南非SAFRICA','恩家那桑NGANASAN','HAZARA-UYGUR-UZBEK','撒丁岛SARDINIA','东西伯利亚EASTSIBERIA','中非CAFRICA','土耳其TURK-伊朗IRAN-高加索CAUCASUS','南亚SSASIA','卡拉什KALASH','PATHAN-SINDHI-BURUSHO','台湾原住民TAIWAN','BALOCHI-MAKRANI-BRAHUI','东亚（核心）EASIA','南美SAMERICA','东北欧NEEUROPE','肯尼亚班图BANTUKENYA','近东NEAREAST','孟德人MENDE','孟加拉BENGALI','北欧NEUROPE','东南亚SEASIA','北非NAFRICA','意大利ITALY','芬兰FINNISH','西南欧SWEUROPE','东亚（日韩）JAPAN-KOREA','俄南TUBALAR','古吉拉特GUJARAT_帕托PATEL','德系犹太人ASHKENAZI','古吉拉特GUJARAT','东非EAFRICA','坦桑尼亚HADZA','中南美CSAMERICA','非洲土著BIAKA','缅甸泰国CAMBODIA-THAI','尼日利亚班图BANTUNIGERIA','塞浦路斯CYPRUS-馬耳他MALTA-西西里島SICILY','冈比亚GAMBIA','中北亚NCASIA','海岛OCEANIA','巴尔干半岛SBALKANS']
#nation_name =  str(nation_name).decode("unicode-escape")
nation_number_total_timed = map(lambda (a, b): a / b, zip(nation_number_total, unweighted_nation_number))
nation_number_total_median =np.median(nation_number_total_timed)
nation_number_total_min_median = [c-nation_number_total_median for c in nation_number_total_timed]
nation_number_total_timed_relu = map(relu,nation_number_total_min_median)
nation_number_total_timed_log = map(psigmoid,nation_number_total_timed_relu)
#nation_number_total_percentage = map(lambda (a, b): a * b, zip(nation_number_total_timed, nation_number_total_min_median))
sum_nub = sum(nation_number_total_timed_log)
result_num = [c/sum_nub for c in nation_number_total_timed_log]
total_list = [nation_name,result_num]

nation_number_total_exped = [n **1000 for n in nation_number_total_timed]
nation_number_total_exped_sum = sum(nation_number_total_exped)
nation_number_total_exped = [n/nation_number_total_exped_sum for n in nation_number_total_exped]
timed_list = [nation_name , nation_number_total_exped]
timed_list = transpose(timed_list)
timed_list_sorted = sorted(timed_list, key = lambda total_list : total_list[1],reverse=True)

print 'Algorithm A 算法A得到的祖源 '
for n in timed_list_sorted:
    for result in n:
        numis = isinstance(result, float)
        if numis is True:
            print('%.4f%%' % ((result) * 100))
        else: print result,
print '----------This is the line 这里是分隔线----------'
print 'Algorithm B 算法B得到的祖源'
total_list = transpose(total_list)
total_list_sorted = sorted(total_list, key = lambda total_list : total_list[1],reverse=True)
for n in total_list_sorted:
    for result in n:
        numis = isinstance(result, float)
        if numis is True:
            print('%.4f%%' % ((result) * 100))
        else: print result,


