#!/usr/bin/env python

from sys import argv, exit
import math


# contributed by Bo Li


if len(argv) != 5:
	print("Usage: python calc_PR.py ntruth ntotal input.ROC output.txt")
	exit(-1)


def output(fout, prog, ntp, nfp, nltp = -1, nlfp = -1):
	""" return delta auc """

	dauc = 0.0
	if nltp < 0:
		recall = 1.0
		precision = ntp * 1.0 / (ntp + nfp)
		fout.write("{}\t{}\t{}\t0\n".format(prog, recall, precision))
	elif ntp == 0 and nfp == 0:
		assert nltp >= 0 and nlfp >= 0 and nltp + nlfp > 0
		lrecall = nltp * 1.0 / ntruth
		lprecision = nltp * 1.0 / (nltp + nlfp)
		recall = 0.0
		precision = lprecision
		if lrecall > 0.0:
			fout.write("{}\t{}\t{}\t0\n".format(prog, recall, precision))
			dauc = 0.5 * lrecall * lprecision
	else:
		recall = ntp * 1.0 / ntruth
		precision = ntp * 1.0 / (ntp + nfp)

		if nltp > ntp:
			lrecall = nltp * 1.0 / ntruth
			lprecision = nltp * 1.0 / (nltp + nlfp)
			
			rate = (nlfp - nfp) * 1.0 / (nltp - ntp)
			trecall = lrecall - 0.01
			x = nltp - ntp - 0.01 * ntruth
			tlrecall = lrecall
			tlprecision = lprecision
			while trecall > recall:
				trecall = (ntp + x) * 1.0 / ntruth
				tprecision = (ntp + x) * 1.0 / (ntp + x + nfp + rate * x)
				fout.write("{}\t{}\t{}\t0\n".format(prog, trecall, tprecision))
				dauc += 0.5 * (tlprecision + tprecision) * (tlrecall - trecall)

				tlrecall = trecall 
				tlprecision = tprecision
				trecall -= 0.01
				x -= 0.01 * ntruth

			dauc += 0.5 * (tlprecision + precision) * (tlrecall - recall)	

		fout.write("{}\t{}\t{}\t1\n".format(prog, recall, precision))
	
	return dauc






ntruth = int(argv[1])
ntotal = int(argv[2])
prog = ""
ltp = lfp = 0
auc = 0.0

with open(argv[3]) as fin, open(argv[4], "w") as fout:
	next(fin)
	for line in fin:
		fields = line.strip().split()

		tp = int(fields[2])
		fp = int(fields[3])

		if prog != fields[0]:
            # prog switch
			if prog != "":
                # process last line of prev prog and report
				auc += output(fout, prog, 0, 0, ltp, lfp)
				print("{}\t{:.2f}".format(prog, auc))
            # first line of next prog, reinit vals
			prog = fields[0]
			ltp = ntruth
			lfp = ntotal - ntruth
			auc = output(fout, prog, ltp, lfp)

        # add to auc
		auc += output(fout, prog, tp, fp, ltp, lfp)
		ltp = tp
		lfp = fp

	if prog != "":
        # last line of file, process last prog results
		auc += output(fout, prog, 0, 0, ltp, lfp)
		print("{}\t{:.2f}".format(prog, auc))
