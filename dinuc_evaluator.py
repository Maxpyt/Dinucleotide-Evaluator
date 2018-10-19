# -*- coding: utf-8 -*-
"""
Created on Wed Feb 07 15:35:20 2018

@author: schnepf
"""

import evaluate_HO
import dinuc_fitter_func as dinuc_fitter
import numpy as np
from scipy.stats import variation
import matplotlib.pyplot as plt
from PyWM_plotter import plotter


# import matplotlib.pyplot as plt

def count_mutations(seq_1, seq_2):
    '''
    function by MvR
    A function that measures the number of mutations between two sequences
    :param seq_1(str): first sequence
    :param seq_2(str): second sequence
    :return (int): number of mutations
    '''
    assert (len(seq_1) == len(seq_2))
    mutation_counter = 0
    for pos in range(len(seq_1)):
        if seq_1[pos] != seq_2[pos]:
            mutation_counter += 1

    return mutation_counter


# a function fetching matching KDs for mutations of certain positions
def calc_0th_order(results, cons, start, stop, KD_col='KD'):
    '''
    function to calculate PWM
    :param results:
    :param cons:
    :param start:
    :param stop:
    :param KD_col:
    :return:
    '''
    prae = cons[:start]
    post = cons[stop:]
    seq = cons[start:stop]
    letters = ['A', 'C', 'G', 'T']
    KD_res = []
    for i in range(len(seq)):
        for j in letters:
            mut = prae + seq[:i] + j + seq[i + 1:] + post
            mut_KD = np.mean(results[results['Sequence'].str.contains(mut)][KD_col])
            KD_res.append(mut_KD)
            # KD_res.append([mut,mut_KD])
    return KD_res


# function to translate KD values to weights (used to generate a PWM)
def KD_to_probs(result_KD, stepsize=16):
    probs = []
    for i in range(len(result_KD) / stepsize):
        sum_rev_prob = np.sum([1 / x for x in result_KD[(i * stepsize):((i + 1) * stepsize)]])
        pos = []
        for j in range(stepsize):
            pos.append((1 / result_KD[i * stepsize + j]) / sum_rev_prob)
        probs.append(pos)
    return probs


# add the predicted affinities based on the linearity assumption
def pred(result_df, cons, start, stop, KD_cons, PWM='', KD_col='KD'):
    result_df = result_df.assign(expected_affinity=0)
    # calculate the PWM based on affinities
    if PWM == '':
        PWM = KD_to_probs(calc_0th_order(result_df, cons, start, stop, KD_col), 4)
    # the PWM_factor is the PWM based on 1 for the consensus (== strongest bound) letter
    PWM_factor = []
    for i in range(len(PWM)):
        PWM_factor.append([x / max(PWM[i]) for x in PWM[i]])
    # iterate over all possible double mutations
    letters = ['A', 'C', 'G', 'T']
    prae = cons[:start]
    post = cons[stop:]
    seq = cons[start:stop]
    # test=[]
    for h in range(len(seq) - 1):
        for i in range(4):
            for j in range(4):
                mut = prae + seq[:h] + letters[i] + letters[j] + seq[h + 2:] + post
                # exp=(1/KD_cons)*PWM_factor[h][i]*PWM_factor[h+1][j]
                # assign the expected affinity based on the linearity assumption  to the dataframe
                # the kd is devided by 1 for the consensus and by a number smaller one for mutations
                exp = KD_cons / PWM_factor[h][i] / PWM_factor[h + 1][j]
                # test.append((mut,letters[i],letters[j],PWM_factor[h][i],PWM_factor[h+1][j]))
                result_df.loc[int(np.where(result_df['Sequence'] == mut)[0]), 'expected_affinity'] = exp  # .values[0]
    return result_df


# define function to create a relative DPWM based on factors assinged in PWM_DPWM() (row by row assingement)
def make_DPWM(result_df, PWM, cons, start, stop):
    letters = ['A', 'C', 'G', 'T']
    prae = cons[:start]
    post = cons[stop:]
    seq = cons[start:stop]
    DPWM = []
    DPWM_abs = []
    for h in range(len(seq) - 1):
        row = []
        kds = []
        for i in range(4):
            for j in range(4):
                # make dinucleotide mutation
                mut = prae + seq[:h] + letters[i] + letters[j] + seq[h + 2:] + post
                row.append(result_df.loc[result_df['Sequence'] == mut, 'factor'].item())
                kds.append(result_df.loc[result_df['Sequence'] == mut, 'KD'].item())
        sum_rev_comp = np.sum([1 / x for x in kds])
        DPWM_abs.append([1 / x / sum_rev_comp for x in kds])
        DPWM.append(row)
    DPWM = np.array(DPWM)
    return DPWM, DPWM_abs


# function to return the PWM and the relative DPWM (factors, deviation from expected)
def PWM_DPWM(result_df, consensus_seq, start, stop):
    PWM = KD_to_probs(calc_0th_order(result_df, consensus_seq, start, stop), 4)
    KD_cons = np.mean(result_df.loc[result_df['Sequence'] == consensus_seq, 'KD'].values)
    # add predicted affinities (lineartiy assumption)
    result_df = pred(result_df, consensus_seq, start, stop, KD_cons, PWM, KD_col='KD')
    # ratio acutual divided by expected
    result_df = result_df.assign(factor=result_df['KD'] / result_df['expected_affinity'])
    DPWM, DPWM_abs = make_DPWM(result_df, PWM, consensus_seq, start, stop)
    return (PWM, DPWM, DPWM_abs)


def find_range(result_df, consensus_seq, col='Sequence'):
    mut_pos = []
    for seq in result_df[col]:
        for i in range(len(seq)):
            if seq[i] != consensus_seq[i]:
                mut_pos.append(i)
            else:
                next
    return (np.min(mut_pos), np.max(mut_pos) + 1)


# function to return the DWPM:
def give_abs_DPWM(PWM, DPWM):
    DPWM_abs = [PWM[0] * 4]
    for i in range(len(PWM[1:])):
        weights = []
        for idx1 in range(4):
            for idx2 in range(4):
                weights.append(PWM[i - 1][idx1] * PWM[i][idx2] * DPWM[i][4 * idx1 + idx2])
        DPWM_abs.append(weights)
    return DPWM_abs


'''
iterations=5
(SSE_out, PWM,DPWM)=dinuc_fitter.main(data_path+dest_name,output_file,iterations,mode = 0,print_results = True)
plt.plot (SSE_out)
'''
###### parameters
data_path = 'P:\\TF-DNA-Binding\\FA\\Data\\180925_Nub_HO_part1_1_MM\\'
data_file = '180925_Nub_HO_part1_1_MM_Results_Kds_clean.txt'
dest_name = data_file[:-3] + 'csv'
# dest_name='171220_Hkb_HO_part1_Kds_augm_mean_based.csv'
filename_order = "P:\\TF-DNA-Binding\\FA\\Oligos\\180607_HO-Nub_part1.xlsx"
pip_scheme = 'P:\TF-DNA-Binding\FA\Oligos\Higher_order-pip-scheme2.csv'
write_csv = True
consensus_seq = 'ACGCCTATGCAAAGGG'
seperator = '\t'
# refs=[10,11,32,33,44,65,66,82]#give reference positions manually -->last position 87- n Water controls
# refs=[10,11,32,33,44,65,66]#last one not complete for Oc
# refs=[10,11,32,33,44,65,66,83]#give reference positions manually -->last position 87- n Water controls for part 2
# refs=[10,11,33,44,65,66,82]# Gt. ref 10 not working
refs = [10, 11, 32, 33, 44, 65, 66, 82]
outliers = []
remove_up_low = False
output_file = data_file[:-4]
set_cons = None
# start=4
# stop=11
random_file = False
name = 'Nub'
output = True

(result_df, KD_ref) = evaluate_HO.main(data_path, data_file, dest_name, filename_order, consensus_seq, write_csv,
                                       pip_scheme, seperator, refs, outliers, remove_up_low=False, return_refs=True)
print KD_ref
print variation(KD_ref)
(start, stop) = find_range(result_df, consensus_seq, col='Sequence')
# assign number of mutations to each sequence (for comparison with older Data, CJ)
result_df.loc[:, 'n_mut'] = [count_mutations(x, consensus_seq) for x in result_df['Sequence']]

if isinstance(set_cons, (int, float)):
    result_df.loc[result_df['Sequence'] == consensus_seq, 'KD'] = set_cons
# construct PWM and DPWM to read them into Marc's motif class

(PWM, DPWM, DPWM_abs) = PWM_DPWM(result_df, consensus_seq, start, stop)
pltr = plotter(PWM, DPWM_abs, DPWM_shift=0)
pltr.plot_PWM()
pltr.plot_DPWM(maxi=2)

motif = dinuc_fitter.Motif(7)
motif.PWM = np.array(PWM) * 4
motif.DPWM = DPWM
# motif.print_probability_file(data_path+output_file+'_on-target')
if output == True:
    motif.print_probability_file(data_path + data_file[:-4] + '_on-target')
    np.savetxt(data_path + output_file + '_on_target.pwm', PWM)
    np.savetxt(data_path + output_file + '_on_target.dpwm', DPWM)
oligomers_all = [(result_df.loc[x, 'Sequence'], result_df.loc[x, 'KD']) for x in range(len(result_df))]
iterations = 20
(SSE_out, PWM, DPWM, iters_best) = dinuc_fitter.main(data_path + dest_name,
                                                     data_path + output_file + '_' + str(iterations) + '_iters',
                                                     iterations, mode=0, print_results=True,
                                                     oligomers_all=oligomers_all, lowest=True, save_results=write_csv)
motif_development = dinuc_fitter.main(data_path + dest_name, data_path + output_file + '_' + str(iterations) + '_iters',
                                      iterations, mode=0, print_results=False, oligomers_all=oligomers_all,
                                      lowest=False, save_results=False)
plt.plot([x[0] for x in motif_development], [x[1] for x in motif_development], '--')
plt.plot([x[0] for x in motif_development], [x[1] for x in motif_development], '*')
plt.xlabel('Iterations')
plt.ylabel('Sum of squared errors')
plt.title(name + ' dinuc_fitter')
pltr = plotter(PWM, DPWM)
pltr.plot_PWM()
pltr.plot_DPWM()
pltr.plot_DPWM()
'''
# possibility to create shuffled matrices for HO Prediction
if isinstance(random_file,str):
    for i in range(100):
            motif.save_random_matrix(random_file+'_'+str(i)+'_first_order_'+name+'.wtmx',name=name,shuffle_positions=True)

'''