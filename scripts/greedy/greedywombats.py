#!/usr/bin/env python3

# sample usage: ./greedyriver.py 512 220 3 0 2
# CSIDH-512 prime
# >=2^220 keys
# B=3
# force the 0 largest primes to be skipped
# try to use 2 cores

from multiprocessing import Pool
import time
import psutil

from fileutils import *

from dacshund import DACsHUND
from dacs32 import all_dacs

from math import log2 , prod
# import scipy.special

from memoized import memoized
import costisog
from optimal_strat import dynamic_programming_algorithm
import sys

from random import randrange

import chain

import json
import csv
import os
from datetime import datetime
from pathlib import Path

first_primes = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, 1381, 1399, 1409, 1423, 1427, 1429]


ells_new_m1 = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743, 751, 757, 761, 769, 773, 787, 797, 809, 811, 821, 823, 827, 829, 839, 853, 857, 859, 863, 877, 881, 883, 887, 907, 911, 919, 929, 937, 941, 947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, 1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, 1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, 1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, 1229, 1231, 1237,1249, 1259, 1277, 1279, 1283,
            1289, 1291, 1297, 1301, 1303, 1307, 1319, 1321, 1327, 1367, 1373, 1381,
            1399, 1409, 1423, 1427, 1429, 1433, 1439,)

cofactors_r = {
    203 : ((4*64+14),[3, 11, 19]),
    192 : ((6*64),[3, 5, 19])}

cofactors_s = {
    209 : ((1*64),[5, 11]),
    195 : ((1*64),[11, 11])}

cofactors_sr = {
    207 : ((1*64),[5, 5, 5]),
    195 : ((1*64),[29])}

cofactors_standard = { 
    226 : ((1*64+4 ),[]),
    225 : ((1*64+4 ),[107]), 
    224 : ((1*64+12),[3,3, 5, 17]),
    223 : ((1*64+26),[3,3, 5]), 
    222 : ((1*64+37),[47]), 
    221 : ((1*64+51),[3]), 
    220 : ((1*64+55),[3,3,3,3, 5]), 
    219 : ((2*64+4 ),[59]), 
    218 : ((2*64+11),[7,7, 13]), 
    217 : ((2*64+24),[3, 5, 7]),
    216 : ((2*64+33),[211]), 
    215 : ((2*64+46),[29]), 
    214 : ((2*64+52),[23, 29]),
    213 : ((3*64),[]),
    212 : ((3*64+18),[]), 
    211 : ((3*64+24),[3, 11]),
    210 : ((3*64+29),[7, 109]),
    209 : ((3*64+43),[73]), 
    208 : ((3*64+53),[73]), 
    207 : ((4*64), []),
    206 : ((4*64+11),[ 3, 13]), 
    205 : ((4*64+19),[ 13, 17]), 
    204 : ((4*64+31),[ 3, 5,5]), 
    203 : ((4*64+39),[ 13, 23]), 
    202 : ((4*64+50),[ 7, 29]), 
    201 : ((5*64), []),
    200 : ((5*64+9 ),[ 3, 13]), 
    199 : ((5*64+16),[ 5, 53]),
    198 : ((5*64+27),[ 191]), 
    197 : ((5*64+35),[ 719]), 
    196 : ((5*64+47),[ 3,3, 31]), 
    195 : ((6*64), []),
    194 : ((6*64+3 ),[ 7,41]), 
    193 : ((6*64+16),[ 67]),
    192 : ((6*64+23),[ 3,3, 41]),
    191 : ((6*64+37),[ 43]),
    190 : ((6*64+47),[ 5, 11]),
    189 : ((6*64+55),[ 5, 31]),
    188 : ((7*64+2 ),[ 7, 19]),
    187 : ((7*64+14),[ 3, 13]),
    186 : ((7*64+23),[ 59]),
    185 : ((7*64+31),[ 7, 43]),
    184 : ((7*64+43),[ 97]),
    183 : ((7*64+55),[ 3, 5]),
    182 : ((8*64), []),
    181 : ((8*64+ 9 ),[ 83]),
    180 : ((8*64+ 19),[ 3,3, 11]),
    179 : ((8*64+ 32),[ 3, 5]),
    178 : ((8*64+ 37),[ 3, 5, 23]),
    177 : ((8*64+ 48),[ 3, 5, 11]),
    176 : ((9*64), []),
    175 : ((9*64+ 3 ),[ 13, 31]),
    174 : ((9*64+ 15),[ 131]),
    173 : ((9*64+ 26),[ 3, 5,5]),
    172 : ((9*64+ 35),[ 107]),
    171 : ((9*64+ 46),[ 47]),
    170 : ((9*64+ 55),[ 137]),
    169 : ((10*64), []),
    168 : ((10*64+ 13),[ 23]),
    167 : ((10*64+ 19),[ 17, 23]),
    166 : ((10*64+ 30),[ 7, 43]),
    165 : ((10*64+ 42),[ 7, 11]),
    164 : ((10*64+ 49),[ 449]),
    163 : ((11*64), []),
    162 : ((11*64+ 9 ),[  3, 11]),
    161 : ((11*64+ 18),[ 41]),
    160 : ((11*64+ 24),[ 3, 7, 41]),
    159 : ((11*64+ 40),[ 13]),
    158 : ((11*64+ 49),[ 23]),
    157 : ((11*64+ 63), []),
    156 : ((12*64+ 3 ),[ 53]),
    155 : ((12*64+ 9 ),[ 5,5, 29]),
    154 : ((12*64+ 21),[ 227]),
    153 : ((12*64+ 31),[ 137]),
    152 : ((12*64+ 41),[ 13,13]),
    151 : ((12*64+ 49),[ 401]),}

cofactors_conf = {
    'standard' : cofactors_standard,
    'reduced'  : cofactors_r,
    'small'    : cofactors_s,
    'small_reduced' : cofactors_sr,}

M = S = 1

def sqrt(p):
  sqrtchain = chain.chain2((p+1)//4)
  sqrtchaincost = chain.cost2(sqrtchain)
  return sqrtchaincost[0]+(sqrtchaincost[1]+1)

def elligator(p):
  return S+M+M+M+S+M+M+sqrt(p)


def inv(p):
  invchain = chain.chain2(p-2)
  invchaincost = chain.cost2(invchain)
  return invchaincost[0]+invchaincost[1]


def dac_search(target,r0,r1,r2,chain,chainlen,best,bestlen):
  if chainlen >= bestlen:
    return best,bestlen
  if r2 > target: return best,bestlen
  if r2<<(bestlen-1-chainlen) < target: return best,bestlen
  if (r2 == target):
    return chain,chainlen
  chain *= 2
  chainlen += 1
  best,bestlen = dac_search(target,r0,r2,r0+r2,chain+1,chainlen,best,bestlen)
  best,bestlen = dac_search(target,r1,r2,r1+r2,chain,chainlen,best,bestlen)
  return best,bestlen

def daclen(target):
  best = None
  bestlen = -1
  while best == None:
    bestlen += 1
    best, bestlen = dac_search(target,1,2,3,0,0,best,bestlen)

  return bestlen

@memoized
def daccost(target):
  return costisog.dac(daclen(target))

def batch_minlenintersec(primes):
    return min(set.intersection(*(DACsHUND[ell] for ell in primes)))

def batch_daclen(primes):
  return batch_minlenintersec(primes)

def batch_daccost(primes):
  return costisog.dac(batch_minlenintersec(primes))

def printstatus(prefix,cost,N0,m0,numprimes1):
  N = N0 if numprimes1 == 0 else N0+(numprimes1,)
  m = m0 if numprimes1 == 0 else m0+(0,)
  print('%s %a %s %s' % (prefix,cost,str(N).replace(' ',''),str(m).replace(' ','')))

@memoized
def batchstart(batchsize):
  B = len(batchsize)
  return [sum(batchsize[:j]) for j in range(B)]


@memoized
def batchstop(batchsize):
  B = len(batchsize)
  return [sum(batchsize[:j+1]) for j in range(B)]


def costfunction(primes0,primes1,N0,m0, strat = False, dummy=False):

  bounds = [sum(N0[:i+1]) for i in range((len(N0)))]
  batch_stop =[b-1 for b in bounds]
  batch_start = [0,]+bounds[:-1]

  keyb = [sum(m0[:i+1]) for i in range(len(m0))]

  key_start = [0,]+keyb[:-1]

  key_stop = [b-1 for b in keyb]

  L_max = []
  L_min = []
  C_mul = []
  for b1, b2, k in zip(batch_start,batch_stop, m0):
      dac = batch_daccost(primes0[b1:b2+1])
      C_mul += [dac]*k
      L_max += [primes0[i] for i in range(b2-k+1, b2+1)]
      L_min += [primes0[i] for i in range(b1, b1+k)]

  C_isog = []
  C_eval = []
  for p_min,  p_max in zip(L_min, L_max):
    ic = 1
    ice = 1
    ic = costisog.isog_matryoshka(p_min, p_max, 0)
    ice = costisog.isog_matryoshka(p_min, p_max, 1) - ic
       

    C_isog += [ic]
    C_eval += [ice]

  C_isog.reverse()
  C_mul.reverse()
  C_eval.reverse()
  L = L_max
  L.reverse()
  cost = dynamic_programming_algorithm(L, C_mul, C_isog, C_eval)

  if strat:
     return cost[1], cost[0]
  return cost[1]


def comb(N, k):
   # from scipy as we don't have it installed
  N = int(N)
  k = int(k)

  if k > N or N < 0 or k < 0:
      return 0

  M = N + 1
  nterms = min(k, N - k)

  numerator = 1
  denominator = 1
  for j in range(1, nterms + 1):
      numerator *= M - j
      denominator *= j

  return numerator // denominator

@memoized
def batchkeys_wombat(x,y):
  ## (-1, 1)
  return comb(x, y)*(2**y)
  ## (-1, 0, +1)
  #return sum([comb(x, i)*(2**i) for i in range(0,y+1)])

def batchkeys_CTIDH(x,y):
  poly = [1]
  for i in range(x):
    newpoly = poly+[0]
    for j in range(len(poly)):
      newpoly[j+1] += poly[j]
    poly = newpoly
  for i in range(y):
    newpoly = poly+[0]
    for j in range(len(poly)):
      newpoly[j+1] += 2*poly[j]
    poly = newpoly
  return poly[x]

@memoized
def keys(N,m):
  result = 1
  for s,b in zip(N,m):
    result *= batchkeys_wombat(s,b)
  return result

# neighboring_intvec; search upwards in non-b directions
def searchdown(minkeyspace,primes0,primes1,N0,m0,cost,b,best):
  if cost >= best[0]:
     return best
  if keys(N0,m0) >= minkeyspace:
    return cost,m0

  return best

def find_initial_batch_sizes(num_batches, primes):
    batch_sizes = [1] * num_batches
    current_batch = 0
    assigned_primes = num_batches
    consecutive_failures = 0
    max_failures = num_batches
    while assigned_primes < len(primes) and consecutive_failures < max_failures:
        test_batch_sizes = batch_sizes[:]
        test_batch_sizes[current_batch] += 1
        if is_dacshund_valid(tuple(test_batch_sizes), num_batches, primes):
            batch_sizes = test_batch_sizes
            assigned_primes += 1
            consecutive_failures = 0
        else:
            consecutive_failures += 1
        current_batch = (current_batch + 1) % num_batches
    if assigned_primes < len(primes):
        return None
    return tuple(batch_sizes)

def is_dacshund_valid(batch_sizes, num_batches, primes):
    prime_index = 0
    for batch_idx in range(num_batches):
        if prime_index >= len(primes):
            break
        batch_size = batch_sizes[batch_idx]
        intersection = set(DACsHUND[primes[prime_index]])
        for prime_offset in range(batch_size):
            current_prime_idx = prime_index + prime_offset
            if current_prime_idx >= len(primes):
                break
            current_prime = primes[current_prime_idx]
            intersection &= DACsHUND[current_prime]
            if not intersection:
                return False
        prime_index += batch_size
    return True

def build_dac_paths(batch_sizes, primes):
    dac_paths = {}
    prime_start_idx = 0
    for batch_idx, batch_size in enumerate(batch_sizes):
        batch_primes = []
        first_prime = primes[prime_start_idx]
        lengths_intersection = DACsHUND[first_prime]
        for offset in range(batch_size):
            prime_idx = prime_start_idx + offset
            if prime_idx < len(primes):
                prime = primes[prime_idx]
                lengths_intersection &= DACsHUND[prime]
                batch_primes.append(prime)
        dac_paths[batch_idx] = {
            "lengths": lengths_intersection,
            "primes": batch_primes
        }
        prime_start_idx += batch_size
    return dac_paths

def steps_by_batchsize(b):
   if b == 1:
      return 1
   if b == 3:
      return 2
   else:
      return b//2

def optimizem(minkeyspace,primes0,primes1,N0,m0=None):
  B0 = len(N0)

  # Idea: randomize this and do multiple runs?
  m0 = [steps_by_batchsize(b) for b in N0]

  # for i, (m, N) in enumerate(zip(m0, N0)):
  #   best = batchkeys_wombat(N, m)
  #   for mnew in range(m, N):q
  #     k =  batchkeys_wombat(N, mnew)
  #     if k > best:
  #       best = k
  #       m0[i] = mnew

  cost = costfunction(primes0,primes1,N0,m0, dyn=True)

  # random runs
  m0_best = []
  cost_best = 1_000_000
  for _ in range(15):
    m0 = [steps_by_batchsize(b)  for b in N0]
    for j in range(4, B0):
      m0[j] -= randrange(5)
    while True:
      best = cost,m0
      for b in range(B0):
        if m0[b] == 0: continue
        newm = list(m0)
        newm[b] -= 1
        newm = tuple(newm)
        newcost = costfunction(primes0,primes1,N0,newm, dyn=True)
        best = searchdown(minkeyspace,primes0,primes1,N0,newm,newcost,b,best)
      if best == (cost,m0): break

      cost, m0 = best
      if cost < cost_best:
        cost_best = cost
        m0_best = m0
        printstatus('improved', cost, N0, tuple(m0), len(primes1))


  return cost_best,m0_best


def optimizeNm(minkeyspace,primes0,primes1,initial_sizes,B,parallelism=1):
  B0 = B #B-1 if len(primes1)>0 else B
  N0 = initial_sizes
  if (N0 == None):
    return None
  #N0 = tuple(len(primes0)//B0+(j<len(primes0)%B0) for j in range(B0))
  cost,m0 = 99999999,[max(1,b//2) for b in N0] #optimizem(minkeyspace,primes0,primes1,N0)
  print(m0)

  while True:
    best = cost,N0,m0
    variants = []
    for b in range(B0):
      if N0[b] <= 1: continue
      for c in range(B0):
        if c == b: continue
        newsize = list(N0)
        newsize[b] -= 1
        newsize[c] += 1
        newsize = tuple(newsize)
        if is_dacshund_valid(newsize,B0,primes0):
          variants += [(minkeyspace,primes0,primes1,newsize,m0)]
    with Pool(parallelism) as p:
      results = p.starmap(optimizem,variants,chunksize=1)
    for (newcost,newm),(_,_,_,newsize,_) in zip(results,variants):
      if newcost < best[0]:
        best = newcost,newsize,newm
    if best == (cost,N0,m0): break
    cost,N0,m0 = best


  return cost,N0,m0


def doit():
  sys.setrecursionlimit(10000)

  cofactor_confs = ['standard','reduced','small','small_reduced']
  cofactor_configuration = 'standard'
  if len(sys.argv) > 1:
    cofactor_configuration = sys.argv[1]
  cofactors = cofactors_conf[cofactor_configuration]

  minkeyspace = 2**221
  # if len(sys.argv) > 2:
  #   minkeyspace = 2**float(sys.argv[2])

  # B = 3
  # if len(sys.argv) > 3:
  #   B = int(sys.argv[3])
  # assert B >= 1
  # assert B <= len(first_primes)

  n_ells = 226
  if len(sys.argv) > 2:
    n_ells = int(sys.argv[2])
  assert 150 <= n_ells <= 226

  parallelism = 24
  if len(sys.argv) > 3:
    parallelism = int(sys.argv[3])

  B_start = 12
  B_end   = 19

  print("start search:  #ell =", n_ells , " with ", parallelism, " threads")
  primes0 = first_primes[:n_ells]
  if cofactor_configuration == 'reduced' or cofactor_configuration == 'small_reduced':
      primes0 = first_primes[1:n_ells]
  primes1 = [first_primes[:1]]

  for B in range(B_start, B_end+1):
    print(f"\n=== Starting optimization for B={B} ===")

    B0 = B
    # Start timing
    start_wall_time = time.time()
    start_cpu_time = get_process_cpu_time()  # Use resource-based measurement

    initial_sizes = find_initial_batch_sizes(B0, primes0)
    print("init N:", initial_sizes)
    Nm = optimizeNm(minkeyspace,primes0,primes1,initial_sizes,B,parallelism)

    # End timing
    end_wall_time = time.time()
    end_cpu_time = get_process_cpu_time()  # Use resource-based measurement

    wall_time = end_wall_time - start_wall_time
    cpu_time = end_cpu_time - start_cpu_time

    if (Nm == None):
      print(f"B={B}: No valid solution found")

    cost,N0,m0 = Nm

    # Print timing information
    print(f"B={B}: Optimization completed")
    # print(f"  Wall time: {format_time(wall_time)}")
    # print(f"  CPU time: {format_time(cpu_time)}")
    # if wall_time > 0:
    #     print(f"  CPU efficiency: {(cpu_time/wall_time)*100:.2f}%")
    # else:
    #     print(f"  CPU efficiency: N/A (wall time too small)")

    # After running optimization
    dac_paths = build_dac_paths(N0, primes0)

    # Save results with timing information
    files = save_optimization_results(
        cost=cost,
        N0=N0,
        m0=m0,
        dac_paths=dac_paths,
        initial_batch_sizes=initial_sizes,
        B=B,
        primes=primes0,
        minkeyspace=minkeyspace,
        cpu_time=cpu_time,
        wall_time=wall_time,
        n_ells=n_ells,
        cofactors_conf=cofactor_configuration
    )
    print()
    print(f"{B}:")
    printstatus(f'[{n_ells}] output',cost,N0,tuple(m0),len(primes1))
    print()

@memoized
def batchstart(batchsize):
  B = len(batchsize)
  return [sum(batchsize[:j]) for j in range(B)]


@memoized
def batchstop(batchsize):
  B = len(batchsize)
  return [sum(batchsize[:j+1]) for j in range(B)]

def wombat_config(N0,m0, primes):
  config = {}
  bounds = [sum(N0[:i+1]) for i in range((len(N0)))]
  config["WOMBATKEYS"] = sum(m0)
  config["batches"] = len(N0)
  config["batch_start"] = [0,]+bounds[:-1]
  config["batch_stop"] = [b-1 for b in bounds]
  keyb = [sum(m0[:i+1]) for i in range(len(m0))]

  config["batch_keybounds_start"] = [0,]+keyb[:-1]

  config["batch_keybounds_stop"] = [b-1 for b in keyb]
  bstart = batchstart(N0)
  bstop = batchstop(N0)

  config["batch_numkeys"] = m0
  config["batch_maxdac"] = [batch_daclen(primes[b1:b2]) for b1,b2 in zip(bstart ,bstop)]
  config["keys"] = log2(keys(N0, m0))
  cost =  costfunction(primes, [], N0, m0, strat=True)
  config["cost"] = cost[0]
  config["strat"] = cost[1]

  dacshund_len = [] 
  for b1,b2 in zip(bstart ,bstop):
     dacshund_len += [batch_daclen(primes[b1:b2])] * len(primes[b1:b2])

  config["primes_dacshund"] = [all_dacs[p][l] for p,l in zip(primes, dacshund_len)]
  return config




if __name__ == '__main__':
  # doit()

  # print("m4l205")
  # skip = 0
  # primes0 = first_primes[skip:194]

  # N = (1, 2, 4, 5, 18, 18, 18, 17, 18, 18, 18, 19, 20, 18)
  # m = (1, 1, 2, 2, 9, 9, 6, 7, 7, 6, 6, 5, 6, 4)
  # config = wombat_config(N, m, primes0)
  # for k,v in config.items():
  #   print(k, ": ", v)

  print("\n\nm4l205 skip 3 old:")
  skip = 1
  primes0 = first_primes[skip:205]
  N = (2, 3, 6, 18, 17, 18, 18, 19, 18, 17, 17, 17, 17, 17)
  m = (1, 2, 3, 9, 8, 8, 6, 7, 7, 4, 4, 4, 4, 3)
  config = wombat_config(N, m, primes0)
  for k,v in config.items():
    print(k, ": ", v)


  print("\n\nm4l205 skip 3 new")
  skip = 1
  primes0 = first_primes[skip:205]
  N = (2, 4, 6, 14, 14, 14, 14, 15, 15, 14, 14, 14, 13, 12, 13, 13, 13)
  m = (1, 2, 3, 7, 7, 6, 6, 7, 6, 6, 4, 5, 4, 3, 2, 2, 2)
  config = wombat_config(N, m, primes0)
  for k,v in config.items():
    print(k, ": ", v)

  # print("\n\nm6l194 skip 3")
  # skip = 1
  # primes0 = first_primes[skip:194]
  # N = (2, 3, 7, 16, 18, 17, 17, 17, 16, 16, 16, 16, 16, 16)
  # m = (1, 2, 3, 8, 9, 8, 8, 6, 7, 5, 4, 5, 4, 3)
  # config = wombat_config(N, m, primes0)
  
  # for k,v in config.items():
  #   print(k, ": ", v)