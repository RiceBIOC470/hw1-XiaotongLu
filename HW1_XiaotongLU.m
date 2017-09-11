% Homework 1. Due before class on 9/5/17

%% Problem 1 - addition with strings

% Fill in the blank space in this section with code that will add 
% the two numbers regardless of variable type. Hint see the matlab
% functions ischar, isnumeric, and str2num. 

%your code should work no matter which of these lines is uncommented. 
x = 3; y = 5; % integers
%x = '3'; y= '5'; %strings
% x = 3; y = '5'; %mixed

%your code goes here
Xiaotong Lu:
for x='3';
y='5';
if ischar(x)
a=str2num(x);
else a=x;
end
if ischar(y)
b=str2num(y);
else b=y;
end
c=a+b;
d=[' x + y = ',num2str(c)];
disp(d)
end
     
 x + y = 8

%output your answer

%% Problem 2 - our first real biology problem. Open reading frames and nested loops.

%part 1: write a piece of code that creates a random DNA sequence of length
% N (i.e. consisting of the letters ATGC) where we will start with N=500 base pairs (b.p.).
% store the output in a variable
% called rand_seq. Hint: the function randi may be useful. 
% Even if you have access to the bioinformatics toolbox, 
% do not use the builtin function randseq for this part. 
Xiaotong Lu:
N = 500; % define sequence length
N=500;
ii=['A','T','C','G'];
num=randi(4,1,N);
rand_seq=ii(num)



%part 2: open reading frames (ORFs) are pieces of DNA that can be
% transcribed and translated. They start with a start codon (ATG) and end with a
% stop codon (TAA, TGA, or TAG). Write a piece of code that finds the longest ORF 
% in your seqeunce rand_seq. Hint: see the function strfind.

Xiaotong Lu:
N=500;
ii=['A','T','C','G'];
num=randi(4,1,N);
rand_seq=ii(num)

rand_seq =

    'AACGGATTCTATAAGTGGACCCCAGATGGTTCGGCTGTGACCATATTCACTAAAGGGCATAATCCCTTCGGCCCACTAAACTCTACTGGTACGGATCCCCCTCAAGTCGAATGTTATTATACACGATCATCAGGGCTCAGCCAATGGCCACGTGATCACCCGGGTTGAATCGCTTGGCTAGCGACTGATTGCATACCAATGCAAGTTGTCCGTCGACTTTGGAGAGGTCTGTACGATCAGGCATCGGGCTGTATCCGGGCGAATCCGCTCCAGATTTCATAATGTTAGTGCTGGTGGTTCGCTCGCCTAAATAGCAAATAACTGTTAAACAGGGTAACGTTCGACCTGTACCTGGGAAACATGCTATGATGACGTACGCGAGTTAGCCGAGTGCGATCGGGCCGAGTACAACTTTTTCAAATGACCGGTTCTGGCTCATTCGGGACAGCGGGTACAATACCCTGCGGTCTCGCAGATCGTCGAAATTGACAGACGTAATT'

k1=strfind(rand_seq,'ATG');
k2=strfind(rand_seq,'TAA');
k3=strfind(rand_seq,'TAG');
k4=strfind(rand_seq,'TGA');
c=[k2,k3,k4];
n=sort(c);
MINUS=bsxfun(@minus,n,k1(:));
Trans=MINUS';
Trans(rem(Trans,3)~=0|Trans<=0)=inf;
ORFMINIMUS=min(Trans);
ORFMinor=ORFMINIMUS(ORFMINIMUS<Inf);
 max(ORFMinor)

%part 3: copy your code in parts 1 and 2 but place it inside a loop that
% runs 1000 times. Use this to determine the probability
% that an sequence of length 500 has an ORF of greater than 50 b.p.


Xiaotong Lu:
for x=1:1000
N=500;
ii=['A','T','C','G'];
num=randi(4,1,N);
rand_seq=ii(num);
k1=strfind(rand_seq,'ATG');
k2=strfind(rand_seq,'TAA');
k3=strfind(rand_seq,'TAG');
k4=strfind(rand_seq,'TGA');
c=[k2,k3,k4];
n=sort(c);
MINUS=bsxfun(@minus,n,k1(:));
Trans=MINUS';
Trans(rem(Trans,3)~=0|Trans<=0)=inf;
ORFMINIMUS=min(Trans);
ORFMinor=ORFMINIMUS(ORFMINIMUS<Inf);
ORFALL=ORFMinor(find(ORFMinor>50));
NUM=int2str(numel(ORFALL));
x=x+1;
NUMALL=int2str(numel(ORFMinor));
end
PROB=sum(NUM(1,:))/1000



%part 4: copy your code from part 3 but put it inside yet another loop,
% this time over the sequence length N. Plot the probability of having an
% ORF > 50 b.p. as a funciton of the sequence length. 
Xiaotong Lu:
NN=300:1000; 
MM=zeros(1,length(NN));
for i=1:length(NN)
for N=300:1000 
for x=1:1000
ii=['A','T','C','G'];
num=randi(4,1,N);
rand_seq=ii(num);
k1=strfind(rand_seq,'ATG');
k2=strfind(rand_seq,'TAA');
k3=strfind(rand_seq,'TAG');
k4=strfind(rand_seq,'TGA');
c=[k2,k3,k4];
n=sort(c);
MINUS=bsxfun(@minus,n,k1(:));
Trans=MINUS';
Trans(rem(Trans,3)~=0|Trans<=0)=inf;
ORFMINIMUS=min(Trans);
ORFMinor=ORFMINIMUS(ORFMINIMUS<Inf);
ORFALL=ORFMinor(find(ORFMinor>50));
NUM=int2str(numel(ORFALL));
x=x+1;
NUMALL=int2str(numel(ORFMinor));
end
PROB=sum(NUM(1,:))/1000;%if compare the ORF greater than50 with the whole ORF, the denominator should be sum(NUMALL(1,:));
MM(i)=PROB;
N=N+1;
end
end
MM;
NN;
plot(NN,MM)
%part 5: Make sure your results from part 4 are sensible. What features
% must this curve have (hint: what should be the value when N is small or when
% N is very large? how should the curve change in between?) Make sure your
% plot looks like this. 
Xiaotong Lu:
The probility will grow gradually from a low number to a greater one
accompanying the increasing of N and when the N is reach to the number great enough,
the probility will be stable and stay unchanged even though the N is continualy increasing.
So the whole curve will look like the letter ¡®S¡¯.
%% problem 3 data input/output and simple analysis

%The file qPCRdata.txt is an actual file that comes from a Roche
%LightCycler qPCR machine. The important columns are the Cp which tells
%you the cycle of amplification and the position which tells you the well
%from the 96 well plate. Each column of the plate has a different gene and
%each row has a different condition. Each gene in done in triplicates so
%columns 1-3 are the same gene, columns 4-6 the same, etc.
%so A1-A3 are gene 1 condition 1, B1-B3 gene 1 condition 2, A4-A6 gene 2
%condition 1, B4-B6 gene2 condition 2 etc. 

% part1: write code to read the Cp data from this file into a vector. You can ignore the last two
% rows with positions beginning with G and H as there were no samples here. 

Xiaotong Lu:
datable=importdata('c:\users\win\desktop\qPCRdata.txt');
data=datable.data;
lat=data(1:72,1)
% Part 2: transform this vector into an array representing the layout of
% the plate. e.g. a 6 row, 12 column array should that data(1,1) = Cp from
% A1, data(1,2) = Cp from A2, data(2,1) = Cp from B1 etc. 
 A=data(1:12,1);
B=data(13:24,1);
C=data(25:36,1);
D=data(37:48,1);
E=data(49:60,1);
F=data(61:72,1);
ALL=[A';B';C';D';E';F']
% Part 3. The 4th gene in columns 10 - 12 is known as a normalization gene.
% That is, it's should not change between conditions and it is used to normalize 
% the expression values for the others. For the other three
% genes, compute their normalized expression in all  conditions, normalized to condition 1. 
% In other words, the fold change between these conditions and condition 1. The
% formula for this is 2^[Cp0 - CpX - (CpN0 - CpNX)] where Cp0 is the Cp for
% the gene in the 1st condition, CpX is the value of Cp in condition X and
% CpN0 and CpNX are the same quantitites for the normalization gene.
% Plot this data in an appropriate way. 
Xiaotong Lu:
DataA=ALL(1:6,1:3);
AverageA=mean(DataA');
DataB=ALL(1:6,4:6);
AverageB=mean(DataB');
DataC=ALL(1:6,7:9);
AverageC=mean(DataC');
DataD=ALL(1:6,10:12);
AverageD=mean(DataD');
Average=[AverageA;AverageB;AverageC;AverageD];
AverageAll=Average';
Nor=zeros(5,3);
for i=1:5
for j=1:3
Nor(i,j)=2^((AverageAll(i+1,j)-AverageAll(1,j))-(AverageAll(i+1,4)-AverageAll(1,4)));
end
end
Uncertainty=zeros(1,3);
Uncertainty=std(Nor)
x=1:5;
NorA=Nor(1:5,1);
NorB=Nor(1:5,2);
NorC=Nor(1:5,3);
meanNor=mean(Nor);
err=zeros(5,3);
for i=1:5
for j=1:3
err(i,j)=Nor(i,j)-meanNor(1,j);
end
end
errorbar(x,NorA,err(1:5,1))
hold on
errorbar(x,NorB,err(1:5,2))
hold on
errorbar(x,NorC,err(1:5,3))



%% Challenge problems that extend the above (optional)

% 1. Write a solution to Problem 2 part 2 that doesn't use any loops at
% all. Hint: start by using the built in function bsxfun to make a matrix of all distances
% between start and stop codons. 

% 2. Problem 2, part 4. Use Matlab to compute the exact solution to this
% problem and compare your answer to what you got previously by testing
% many sequences. Plot both on the same set of axes. Hint: to get started 
% think about the following:
% A. How many sequences of length N are there?
% B. How many ways of making an ORF of length N_ORF are there?
% C. For each N_ORF how many ways of position this reading frame in a
% sequence of length N are there?

% 3. Problem 3. Assume that the error in each Cp is the standard deviation
% of the three measurements. Add a section to your code that propogates this
% uncertainty to the final results. Add error bars to your plot. (on
% propagation of error, see, for example:
% https://en.wikipedia.org/wiki/Propagation_of_uncertainty


