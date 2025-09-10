#!/usr/bin/env python3

M = 1
S = 1

from memoized import memoized
from costpoly import tree1,multiprod2,multiprod2_selfreciprocal,multieval_precompute,multieval_postcompute

x2 = S+S
x2DBL = M+M+M+M # but save 1M if affine
xDBL = x2+x2DBL
xADD = M+M+S+S+M+M

def dac(daclen):
  return xDBL+xADD+daclen*xADD

@memoized
def isog(lmax,push,bsgs, lmin = -1):
  bs,gs = bsgs

  assert lmax > 0
  #assert lmax%2
  lbits = 6
  while lmax>>lbits: lbits += 2

  result = 0

  br = (lmax-1) // 2
  bl = (lmin-1) // 2
  t = br - bl


  assert bs >= 0
  assert gs >= 0

  if bs == 0 or gs == 0:
    result += xDBL
    result += (br - 1)*xADD
    result += (br)*2*M       # A
    result += (br)*4*M*push  # push points
    if lmin != -1:
      for i in range(bl+1, br+1):
        for _ in range(bl+1-t, i):
          result += 2*M # proj eq overhead for dummy freeness
  
  else: # velusqrt
    multiples = set()
    multiples.add(2) # ???
    for j in range(3,2*bs+1,2): multiples.add(j)  # 490 - 507
    for j in range(4,lmax+1-4*bs*gs,2): multiples.add(j) # 512 - 529
    multiples.add(2*bs)  #
    multiples.add(4*bs)  # ??
    for j in range(6*bs,2*bs*(2*gs+1),4*bs): multiples.add(j) # 523 - 572 ??
  
    result += len(multiples)*xADD # actually xADD or xDBL

    result += tree1[bs]                                 # 623
    result += gs*(3*S+4*M) # biquad_curve               # 632
    result += push*(3*S+M) # biquad_precompute_point
    result += push*gs*6*M # biquad_postcompute_point
    result += push*multiprod2[gs] # pushing point
    result += 2*multiprod2_selfreciprocal[gs] # pushing curve 
    result += multieval_precompute[bs][gs*2+1] # reciprocal of root of product tree
    result += 2*(push+1)*multieval_postcompute[bs][gs*2+1] # scaled remainder tree
    result += 2*M*(push+1)*(bs-1) # accumulating multieval results
    result += (2*M+(4*M*push))*(br-2*bs*gs) # stray points at the end

    result += t*(xDBL + 2*M) # check xDBLL([i]K) == [2i]K
    for _ in range(bs, ((br-2*bs*gs) // 2)):
      result += xADD #compute missing points
  
  result += push*(S+S+M+M) # final point evaluation
  result += 2*(S+M+M+(S+S+M)*(lbits//2-1)) # powpow8mod
  return result


def optimize_bs_gs(l):
  bsdf = {3: (0, 0), 5: (0, 0), 7: (0, 0), 11: (0, 0), 13: (0, 0), 17: (0, 0), 19: (0, 0), 23: (0, 0), 29: (0, 0), 31: (0, 0), 37: (0, 0), 41: (0, 0), 43: (0, 0), 47: (0, 0), 53: (0, 0), 59: (0, 0), 61: (0, 0), 67: (0, 0), 71: (0, 0), 73: (0, 0), 79: (0, 0), 83: (0, 0), 89: (0, 0), 97: (6, 4), 101: (6, 4), 103: (6, 4), 107: (6, 4), 109: (6, 4), 113: (6, 4), 127: (6, 5), 131: (8, 4), 137: (8, 4), 139: (8, 4), 149: (6, 6), 151: (6, 6), 157: (6, 6), 163: (8, 5), 167: (8, 5), 173: (8, 5), 179: (8, 5), 181: (8, 5), 191: (8, 5), 193: (8, 6), 197: (8, 6), 199: (8, 6), 211: (8, 6), 223: (8, 6), 227: (8, 7), 229: (8, 7), 233: (8, 7), 239: (8, 7), 241: (10, 6), 251: (10, 6), 257: (8, 8), 263: (8, 8), 269: (8, 8), 271: (8, 8), 277: (8, 8), 281: (10, 7), 283: (10, 7), 293: (12, 6), 307: (12, 6), 311: (12, 6), 313: (12, 6), 317: (12, 6), 331: (10, 8), 337: (14, 6), 347: (14, 6), 349: (14, 6), 353: (14, 6), 359: (14, 6), 367: (10, 9), 373: (10, 9), 379: (10, 9), 383: (10, 9), 389: (12, 8), 397: (14, 7), 401: (14, 7), 409: (14, 7), 419: (14, 7), 421: (14, 7), 431: (14, 7), 433: (12, 9), 439: (12, 9), 443: (12, 9), 449: (14, 8), 457: (14, 8), 461: (14, 8), 463: (14, 8), 467: (14, 8), 479: (14, 8), 487: (12, 10), 491: (12, 10), 499: (12, 10), 503: (12, 10), 509: (14, 9), 521: (14, 9), 523: (14, 9), 541: (14, 9), 547: (14, 9), 557: (14, 9), 563: (14, 10), 569: (14, 10), 571: (14, 10), 577: (14, 10), 587: (14, 10), 593: (14, 10), 599: (14, 10), 601: (14, 10), 607: (14, 10), 613: (14, 10), 617: (14, 11), 619: (14, 11), 631: (14, 11), 641: (16, 10), 643: (16, 10), 647: (16, 10), 653: (16, 10), 659: (16, 10), 661: (16, 10), 673: (14, 12), 677: (14, 12), 683: (14, 12), 691: (14, 12), 701: (14, 12), 709: (16, 11), 719: (16, 11), 727: (16, 11), 733: (14, 13), 739: (14, 13), 743: (14, 13), 751: (14, 13), 757: (14, 13), 761: (14, 13), 769: (16, 12), 773: (16, 12), 787: (16, 12), 797: (16, 12), 809: (16, 12), 811: (16, 12), 821: (16, 12), 823: (16, 12), 827: (16, 12), 829: (16, 12), 839: (16, 13), 853: (16, 13), 857: (16, 13), 859: (16, 13), 863: (16, 13), 877: (18, 12), 881: (18, 12), 883: (18, 12), 887: (18, 12), 907: (16, 14), 911: (16, 14), 919: (16, 14), 929: (16, 14), 937: (18, 13), 941: (18, 13), 947: (18, 13), 953: (18, 13), 967: (20, 12), 971: (20, 12), 977: (20, 12), 983: (20, 12), 991: (20, 12), 997: (20, 12), 1009: (18, 14), 1013: (18, 14), 1019: (18, 14), 1021: (18, 14), 1031: (16, 16), 1033: (16, 16), 1039: (16, 16), 1049: (20, 13), 1051: (20, 13), 1061: (22, 12), 1063: (22, 12), 1069: (22, 12), 1087: (18, 15), 1091: (18, 15), 1093: (18, 15), 1097: (18, 15), 1103: (18, 15), 1109: (18, 15), 1117: (18, 15), 1123: (20, 14), 1129: (20, 14), 1151: (22, 13), 1153: (22, 13), 1163: (18, 16), 1171: (18, 16), 1181: (18, 16), 1187: (18, 16), 1193: (18, 16), 1201: (20, 15), 1213: (20, 15), 1217: (20, 15), 1223: (20, 15), 1229: (18, 17), 1231: (18, 17), 1237: (22, 14), 1249: (22, 14), 1259: (22, 14), 1277: (22, 14), 1279: (22, 14), 1283: (20, 16), 1289: (20, 16), 1291: (20, 16), 1297: (20, 16), 1301: (20, 16), 1303: (20, 16), 1307: (20, 16), 1319: (20, 16), 1321: (22, 15), 1327: (22, 15), 1367: (20, 17), 1373: (20, 17), 1381: (20, 17), 1399: (20, 17), 1409: (22, 16), 1423: (22, 16), 1427: (22, 16), 1429: (22, 16), 1433: (22, 16), 1439: (22, 16)}
  
  return bsdf[l]

@memoized
def optimize(l,push):
  #bsgs = optimize_bs_gs(l)
  #return isog(l,push, bsgs), bsgs
  
  #tmp = bsdf[l]

  best,bestbsgs = isog(l,push,(0,0)),(0,0)

  # XXX: precompute more tree1 etc.; extend bs,gs limits
  for bs in range(2,33,2):
    for gs in range(1,2*bs+1):
      if 2*bs*gs > (l-1)//2: break
      if gs >= 32: break
      if bs > 3*gs: continue
      result = isog(l,push,(bs,gs))
      if result < best:
        best,bestbsgs = result,(bs,gs)

  return best,bestbsgs


def isog_matryoshka(l_min, l_max, push):
  (bs, gs) = optimize_bs_gs(l_min)

  velu = isog(l_max, push, (bs,gs), l_min)
  return velu
