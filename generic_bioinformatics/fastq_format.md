# Fastq Files Format

Each entry in a FASTQ file consists of four lines:
- Sequence identifier
	- ```@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos> <read>:<is filtered>:<control number>:<sample number>```
- Sequence
- Quality score identifier line (consisting only of a +)
- Quality score

*For undetermined FASTQ*
In the sequence identifier the sequence observed in the index read is written to the FASTQ header in place of the sample number


Example entry:
```
@SIM:1:FCX:1:15:6329:1045 1:N:0:2
TCGCACTCAACGCCCTGCATATGACAAGACAGAATC
+
<>;##=><9=AAAAAAAAAA9#:<#<;<<<????#=
```

Sequence identifier entries:
- @
	+ Each sequence identifier line starts with @
- <instrument>
	+ Instrument id ([a-zA-Z0-9])
- <run number>
	+ Run number on instrument ([0-9])
- <flowcell ID> 
	+ Flowcell ID ([a-zA-Z0-9])
- <lane>
	+ Lane number
- <tile>
	+ Tile number 
- <x_pos>
	+ X coordinate of cluster
- <y_pos>
	+ Y coordinate of cluster
- <read>
	+ Read number. 1 is single read, 2 is paired-end
- <is filtered>
	+ Y if it is filtered (did not pass), N otherwise
- <control number> 
	+ 0 when none of the control bits are on, otherwise it is an even number. On HiSeq X systems, control specification is not performed and this number is always 0.
- <sample number>
	+ Sample number from sample sheet
	

# Quality scores
Express an error probability
Q(A) = -10 Log10(P(~A))
where P(~A) is the estimated probability of an assertion A being wrong

Quality score Q(A) - Error probability P(~A)
10 -> 0.1
20 -> 0.01
30 -> 0.001

## Quality score encoding
ASCII code +33
```
Symbol	ASCII Code	Q-Score
!	33	0
"	34	1
#	35	2
$	36	3
%	37	4
&	38	5
'	39	6
(	40	7
)	41	8
*	42	9
+	43	10
,	44	11
-	45	12
.	46	13
/	47	14
0	48	15
1	49	16
2	50	17
3	51	18
4	52	19
5	53	20
6	54	21
7	55	22
8	56	23
9	57	24
:	58	25
;	59	26
<	60	27
=	61	28
>	62	29
?	63	30
@	64	31
A	65	32
B	66	33
C	67	34
D	68	35
E	69	36
F	70	37
G	71	38
H	72	39
I	73	40
```










