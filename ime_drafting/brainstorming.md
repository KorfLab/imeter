**Brainstorming!**

*Different ways to improve imeter:*

+ account for degenerative sequences
+ account for – and + strands
+ set up a way to weigh the beginning of the intron vs the end of the intron (should work on finding a geometric decay in it)
	+ would be calculated during counts, so like count += decay, total += decay
+ account for sequence repetition due to splice variants
+ decoding, putting things into libraries (proper program in argparse, put functions into some shared library)
	+ trainer program
	+ decoder program
+ scientific studies
+ controls? Experimental data and controls ( R squared value
+ do certain kmers predict long introns? Long introns will contribute more kmers, if one is more common in long kmers, will be disproportionately weighted.
	+ could weight kmers on length of parent intron
	+ could weight kmers based on composition of sequence? might be something good
+ use expression as a response variable,  like a checker.
+ incorporate 4-5 fold cross val, eventually*
+ sort sequences at end, by score, - to +
+ we could score the kmers separately??
	+ if we score without weight based on position in intron?
	+ then after with weight? or developing a secondary metric?
		+ like adding a second log-odds, across the single intron cutoff at an arbitrary(?) point, sortof the same math ? ? ?
	+ is there a virtue in more stringent cutoffs? i.e, is there a virtue in making a deadzone, so to speak, where introns are ignored because they aren't 'distal' or 'proximal' enough?
	+ in general, the approach relies on big data sets, is there something we could account for for the future for smaller ones (this is out of scope for a rotation, more out of curiosity)
	+ split based on expression level?
		+ could expression be a separate function? perhaps accounting for constitutive expression or high expression
		+ what is the optimal way to split? there are many parameters
+ potential to experimentally determine best split options....
+ investigate python itertools combinations, or something to handle permutations, ie (ACTG, 5), then it makes all the possible combinations
+ separate training and testing into two different programs, allowing the testing program to use an established "training file"


can't be flat...can it? every 50 bases you turn it down by like 0.1? because
depending on how long or short the intron is you won't account for the start vs. end consistently
so maybe we take the whole length and every increment we decrement the decay. This makes it
consistent across introns regardless of length.

In order to do this, I think I need a cutoff, like where in the intron the significance drops
off. May have to be experimentally determined...but if it's a cutoff, a flat value may make more
sense. Does a 1000bp intron have a larger "range of significance" than a 100bp one? Good to ask!

decay function draft, flat decrement value
lets assume without cause that anything after 100 bp is less important...by some...proportion...
sig_assumed = 100
rate_assumed = 0.5
if i+k <= sig_assumed -1:
		decay *= rate_assumed
		sig_assumed += sig_assumed
okay...now lets say instead the decay is more based on how far on any given intron you've gone down.
The assumptions here are again a suitable rate of decay and how often you want to decrement.
rate_assumed = 0.5
proportion_assumed = 0.1
sig_interval = len(s) * proportion_assumed
#then again, the same idea
if i+k <= sig_interval -1:
		decay *= rate_assumed
		sig_interval += sig_interval


		AFTER 1000BP, IT DOES A LOT LESS! SO perhaps set a high limit,
		then make sure that decay makes kmers count for little to nothing at that point.

		output like a standard table. kmer, and their frequencies
