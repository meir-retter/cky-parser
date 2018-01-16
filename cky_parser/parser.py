#!/usr/bin/python3

import json
import os, sys


def make_dics(counts_file_name):
	wordcounts = {} # note that rare words themselves will not appear here
	ntcounts = {} # nonterminal counts
	bncounts = {} # binary rule counts. bncounts[(X, Y1, Y2)] = num times X -> Y1Y2 is used
	uncounts = {} # unary rule counts. uncounts([(X, w)]) = num times X -> w was used
	counts_file = open(counts_file_name, "r")
	for line in [line.split() for line in counts_file.readlines()]:
		if line[1] == "UNARYRULE":
			word = line[3]
			num = int(line[0])
			if word in wordcounts.keys():
				wordcounts[word] += num
			else:
				wordcounts[word] = num
			# also,
			first = line[2]
			uncounts[(first, word)] = num
		elif line[1] == "NONTERMINAL":
			nt = line[2]
			num = int(line[0])
			ntcounts[nt] = num
		elif line[1] == "BINARYRULE":
			first = line[2]
			second = line[3]
			third = line[4]
			num = int(line[0])
			bncounts[(first, second, third)] = num
	is_rare = {}
	for word in wordcounts.keys():
		if wordcounts[word] < 5:
			is_rare[word] = True
		else:
			is_rare[word] = False
	return wordcounts, ntcounts, bncounts, uncounts, is_rare


def q_bn(X, Y1, Y2, ntcounts, bncounts):
	return bncounts[(X, Y1, Y2)] / ntcounts[X]

def q_un(X, w, wordcounts, ntcounts, uncounts):
	if w in wordcounts.keys():
		return uncounts[(X, w)] / ntcounts[X]
	else:
		return uncounts[(X, "_RARE_")] / ntcounts[X]

def CKY(sentence, wordcounts, ntcounts, bncounts, uncounts):
	#takes input as list of strings. we insert a dummy var so indexing begins at 1
	memo = {} # this is pi, the dictionary
	bp = {}
	x = ["DUMMY"] + sentence # this way, x[1] returns x1
	n = len(x) - 1
	nonterminals = ntcounts.keys()
	bnrules = bncounts.keys()
	for i in range(1,n+1):
		for X in nonterminals:
			if ((X, x[i]) in uncounts.keys()) or ((x[i] not in wordcounts.keys()) and (X, "_RARE_") in uncounts.keys()):
				memo[(i, i, X)] = q_un(X, x[i], wordcounts, ntcounts, uncounts)
			else:
				memo[(i, i, X)] = 0
	for l in range(1,n):
		for i in range(1,n-l+1):
			j = i+l
			for X in nonterminals:
				max_pi_val = 0
				best_bin_rule = None
				best_split_point = None
				Xrules = [rule for rule in bnrules if rule[0] == X]
				for rule in Xrules:
					Y = rule[1]
					Z = rule[2]
					for split_point in range(i, j):
						if (i, split_point, Y) not in memo.keys() or (split_point+1, j, Z) not in memo.keys():
							pi_val = 0
						else:
							pi_val = q_bn(X, Y, Z, ntcounts, bncounts) * memo[(i, split_point, Y)] * memo[(split_point+1, j, Z)]
						if pi_val > max_pi_val:
							max_pi_val = pi_val
							best_bin_rule = (X, Y, Z)
							best_split_point = split_point
				if max_pi_val > 0:
					memo[(i, j, X)] = max_pi_val
					bp[(i, j, X)] = (best_bin_rule, best_split_point)
	if (1, n, "S") in memo.keys() and memo[(1, n, "S")] != 0:
		return memo[(1, n, "S")], bp, "S"
	else:
		max_pi_val = 0
		bestX = None
		for X in nonterminals:
			if (1, n, X) not in memo.keys():
				pi_val = 0
			else:
				pi_val = memo[(1, n, X)]
			if pi_val > max_pi_val:
				max_pi_val = pi_val
				bestX = X
		if max_pi_val > 0:
			return max_pi_val, bp, bestX
		else:
			return 0, None, None

def print_sentence(sentence):
	for word in sentence:
		print(word, end=' ')




def make_tree(sentence, i, j, bp, root_nt="S"):
	# now will not use dummy variable, so
	# will need to take into account 0-indexing
	if i == j:
		return [root_nt, sentence[i-1]]
	else:
		best_bin_rule, best_split_point = bp[(i, j, root_nt)]
		left_subtree = make_tree(sentence, i, best_split_point, bp, best_bin_rule[1])
		right_subtree = make_tree(sentence, best_split_point+1, j, bp, best_bin_rule[2])
		return [root_nt, left_subtree, right_subtree]


def make_rare(tree, is_rare):
	if len(tree) == 2:
		if type(tree[1]) is str:
			if is_rare[tree[1]]:
				tree[1] = "_RARE_"
		else:
			make_rare(tree[1], is_rare)
	elif len(tree) == 3:
		if type(tree[1]) is str:
			if is_rare[tree[1]]:
				tree[1] = "_RARE_"
		else:
			make_rare(tree[1], is_rare)
		if type(tree[2]) is str:
			if is_rare[tree[2]]:
				tree[2] = "_RARE_"
		else:
			make_rare(tree[2], is_rare)


def send_training_data_to_output_file(counts_file_name, training_file_name, output_file_name):
	traindat = open(training_file_name, "r")
	is_rare = make_dics(counts_file_name)[4]
	trees = []
	for line in traindat:
		trees.append(json.loads(line))

	for tree in trees:
		make_rare(tree, is_rare)
	parse_train_with_rare = open(output_file_name, "w")
	for tree in trees:
		parse_train_with_rare.write(json.dumps(tree)+'\n')



def do_parsing_and_write_trees_to_pred_file(counts_file_name, output_file_name):
	wordcounts, ntcounts, bncounts, uncounts, is_rare = make_dics(counts_file_name)
	parse_dev = open("parse_dev.dat", "r")
	pred_file = open(output_file_name, "w")
	sentences = [line.split() for line in parse_dev.readlines()]
	for sentence in sentences:
		val, bp, root_nt = CKY(sentence, wordcounts, ntcounts, bncounts, uncounts)
		pred_file.write(json.dumps(make_tree(sentence, 1, len(sentence), bp, root_nt))+'\n')

if __name__ == "__main__":
	training_file_name = sys.argv[2]
	q_num = sys.argv[1]
	os.system("python count_cfg_freq.py " + training_file_name + " > cfg.counts")
	counts_file_name = "cfg.counts"
	if q_num == "q4":
		output_file_name = sys.argv[3]
		send_training_data_to_output_file(counts_file_name, training_file_name, output_file_name)

	else:
		output_file_name = sys.argv[4]
		do_parsing_and_write_trees_to_pred_file(counts_file_name, output_file_name)
