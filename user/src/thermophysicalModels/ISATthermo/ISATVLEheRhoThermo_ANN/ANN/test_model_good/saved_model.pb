��
��
^
AssignVariableOp
resource
value"dtype"
dtypetype"
validate_shapebool( �
~
BiasAdd

value"T	
bias"T
output"T" 
Ttype:
2	"-
data_formatstringNHWC:
NHWCNCHW
8
Const
output"dtype"
valuetensor"
dtypetype
.
Identity

input"T
output"T"	
Ttype
q
MatMul
a"T
b"T
product"T"
transpose_abool( "
transpose_bbool( "
Ttype:

2	
>
Maximum
x"T
y"T
z"T"
Ttype:
2	
�
MergeV2Checkpoints
checkpoint_prefixes
destination_prefix"
delete_old_dirsbool("
allow_missing_filesbool( �

NoOp
M
Pack
values"T*N
output"T"
Nint(0"	
Ttype"
axisint 
C
Placeholder
output"dtype"
dtypetype"
shapeshape:
@
ReadVariableOp
resource
value"dtype"
dtypetype�
@
RealDiv
x"T
y"T
z"T"
Ttype:
2	
E
Relu
features"T
activations"T"
Ttype:
2	
o
	RestoreV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0�
l
SaveV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0�
?
Select
	condition

t"T
e"T
output"T"	
Ttype
H
ShardedFilename
basename	
shard

num_shards
filename
-
Sqrt
x"T
y"T"
Ttype:

2
�
StatefulPartitionedCall
args2Tin
output2Tout"
Tin
list(type)("
Tout
list(type)("	
ffunc"
configstring "
config_protostring "
executor_typestring ��
@
StaticRegexFullMatch	
input

output
"
patternstring
N

StringJoin
inputs*N

output"
Nint(0"
	separatorstring 
<
Sub
x"T
y"T
z"T"
Ttype:
2	
�
VarHandleOp
resource"
	containerstring "
shared_namestring "
dtypetype"
shapeshape"#
allowed_deviceslist(string)
 �"serve*2.11.02v2.11.0-rc2-17-gd5b57ca93e58��
f
ConstConst*
_output_shapes

:*
dtype0*)
value B"��aL���4�r<�r<
h
Const_1Const*
_output_shapes

:*
dtype0*)
value B"
�UF�:;L?�O>
^
countVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_namecount
W
count/Read/ReadVariableOpReadVariableOpcount*
_output_shapes
: *
dtype0
^
totalVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nametotal
W
total/Read/ReadVariableOpReadVariableOptotal*
_output_shapes
: *
dtype0
b
count_1VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	count_1
[
count_1/Read/ReadVariableOpReadVariableOpcount_1*
_output_shapes
: *
dtype0
b
total_1VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	total_1
[
total_1/Read/ReadVariableOpReadVariableOptotal_1*
_output_shapes
: *
dtype0
b
count_2VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	count_2
[
count_2/Read/ReadVariableOpReadVariableOpcount_2*
_output_shapes
: *
dtype0
b
total_2VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	total_2
[
total_2/Read/ReadVariableOpReadVariableOptotal_2*
_output_shapes
: *
dtype0
b
count_3VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	count_3
[
count_3/Read/ReadVariableOpReadVariableOpcount_3*
_output_shapes
: *
dtype0
b
total_3VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	total_3
[
total_3/Read/ReadVariableOpReadVariableOptotal_3*
_output_shapes
: *
dtype0
b
count_4VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	count_4
[
count_4/Read/ReadVariableOpReadVariableOpcount_4*
_output_shapes
: *
dtype0
b
total_4VarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	total_4
[
total_4/Read/ReadVariableOpReadVariableOptotal_4*
_output_shapes
: *
dtype0
r
Adam/v/c/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_nameAdam/v/c/bias
k
!Adam/v/c/bias/Read/ReadVariableOpReadVariableOpAdam/v/c/bias*
_output_shapes
:*
dtype0
r
Adam/m/c/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_nameAdam/m/c/bias
k
!Adam/m/c/bias/Read/ReadVariableOpReadVariableOpAdam/m/c/bias*
_output_shapes
:*
dtype0
z
Adam/v/c/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
: * 
shared_nameAdam/v/c/kernel
s
#Adam/v/c/kernel/Read/ReadVariableOpReadVariableOpAdam/v/c/kernel*
_output_shapes

: *
dtype0
z
Adam/m/c/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
: * 
shared_nameAdam/m/c/kernel
s
#Adam/m/c/kernel/Read/ReadVariableOpReadVariableOpAdam/m/c/kernel*
_output_shapes

: *
dtype0
v
Adam/v/phi/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:* 
shared_nameAdam/v/phi/bias
o
#Adam/v/phi/bias/Read/ReadVariableOpReadVariableOpAdam/v/phi/bias*
_output_shapes
:*
dtype0
v
Adam/m/phi/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:* 
shared_nameAdam/m/phi/bias
o
#Adam/m/phi/bias/Read/ReadVariableOpReadVariableOpAdam/m/phi/bias*
_output_shapes
:*
dtype0
~
Adam/v/phi/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
: *"
shared_nameAdam/v/phi/kernel
w
%Adam/v/phi/kernel/Read/ReadVariableOpReadVariableOpAdam/v/phi/kernel*
_output_shapes

: *
dtype0
~
Adam/m/phi/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
: *"
shared_nameAdam/m/phi/kernel
w
%Adam/m/phi/kernel/Read/ReadVariableOpReadVariableOpAdam/m/phi/kernel*
_output_shapes

: *
dtype0
r
Adam/v/P/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_nameAdam/v/P/bias
k
!Adam/v/P/bias/Read/ReadVariableOpReadVariableOpAdam/v/P/bias*
_output_shapes
:*
dtype0
r
Adam/m/P/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_nameAdam/m/P/bias
k
!Adam/m/P/bias/Read/ReadVariableOpReadVariableOpAdam/m/P/bias*
_output_shapes
:*
dtype0
z
Adam/v/P/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
: * 
shared_nameAdam/v/P/kernel
s
#Adam/v/P/kernel/Read/ReadVariableOpReadVariableOpAdam/v/P/kernel*
_output_shapes

: *
dtype0
z
Adam/m/P/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
: * 
shared_nameAdam/m/P/kernel
s
#Adam/m/P/kernel/Read/ReadVariableOpReadVariableOpAdam/m/P/kernel*
_output_shapes

: *
dtype0
r
Adam/v/T/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_nameAdam/v/T/bias
k
!Adam/v/T/bias/Read/ReadVariableOpReadVariableOpAdam/v/T/bias*
_output_shapes
:*
dtype0
r
Adam/m/T/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_nameAdam/m/T/bias
k
!Adam/m/T/bias/Read/ReadVariableOpReadVariableOpAdam/m/T/bias*
_output_shapes
:*
dtype0
z
Adam/v/T/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
: * 
shared_nameAdam/v/T/kernel
s
#Adam/v/T/kernel/Read/ReadVariableOpReadVariableOpAdam/v/T/kernel*
_output_shapes

: *
dtype0
z
Adam/m/T/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
: * 
shared_nameAdam/m/T/kernel
s
#Adam/m/T/kernel/Read/ReadVariableOpReadVariableOpAdam/m/T/kernel*
_output_shapes

: *
dtype0
t
Adam/v/12/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nameAdam/v/12/bias
m
"Adam/v/12/bias/Read/ReadVariableOpReadVariableOpAdam/v/12/bias*
_output_shapes
: *
dtype0
t
Adam/m/12/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nameAdam/m/12/bias
m
"Adam/m/12/bias/Read/ReadVariableOpReadVariableOpAdam/m/12/bias*
_output_shapes
: *
dtype0
|
Adam/v/12/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@ *!
shared_nameAdam/v/12/kernel
u
$Adam/v/12/kernel/Read/ReadVariableOpReadVariableOpAdam/v/12/kernel*
_output_shapes

:@ *
dtype0
|
Adam/m/12/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@ *!
shared_nameAdam/m/12/kernel
u
$Adam/m/12/kernel/Read/ReadVariableOpReadVariableOpAdam/m/12/kernel*
_output_shapes

:@ *
dtype0
r
Adam/v/9/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nameAdam/v/9/bias
k
!Adam/v/9/bias/Read/ReadVariableOpReadVariableOpAdam/v/9/bias*
_output_shapes
: *
dtype0
r
Adam/m/9/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nameAdam/m/9/bias
k
!Adam/m/9/bias/Read/ReadVariableOpReadVariableOpAdam/m/9/bias*
_output_shapes
: *
dtype0
z
Adam/v/9/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@ * 
shared_nameAdam/v/9/kernel
s
#Adam/v/9/kernel/Read/ReadVariableOpReadVariableOpAdam/v/9/kernel*
_output_shapes

:@ *
dtype0
z
Adam/m/9/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@ * 
shared_nameAdam/m/9/kernel
s
#Adam/m/9/kernel/Read/ReadVariableOpReadVariableOpAdam/m/9/kernel*
_output_shapes

:@ *
dtype0
r
Adam/v/6/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nameAdam/v/6/bias
k
!Adam/v/6/bias/Read/ReadVariableOpReadVariableOpAdam/v/6/bias*
_output_shapes
: *
dtype0
r
Adam/m/6/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nameAdam/m/6/bias
k
!Adam/m/6/bias/Read/ReadVariableOpReadVariableOpAdam/m/6/bias*
_output_shapes
: *
dtype0
z
Adam/v/6/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@ * 
shared_nameAdam/v/6/kernel
s
#Adam/v/6/kernel/Read/ReadVariableOpReadVariableOpAdam/v/6/kernel*
_output_shapes

:@ *
dtype0
z
Adam/m/6/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@ * 
shared_nameAdam/m/6/kernel
s
#Adam/m/6/kernel/Read/ReadVariableOpReadVariableOpAdam/m/6/kernel*
_output_shapes

:@ *
dtype0
r
Adam/v/3/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nameAdam/v/3/bias
k
!Adam/v/3/bias/Read/ReadVariableOpReadVariableOpAdam/v/3/bias*
_output_shapes
: *
dtype0
r
Adam/m/3/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_nameAdam/m/3/bias
k
!Adam/m/3/bias/Read/ReadVariableOpReadVariableOpAdam/m/3/bias*
_output_shapes
: *
dtype0
z
Adam/v/3/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@ * 
shared_nameAdam/v/3/kernel
s
#Adam/v/3/kernel/Read/ReadVariableOpReadVariableOpAdam/v/3/kernel*
_output_shapes

:@ *
dtype0
z
Adam/m/3/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@ * 
shared_nameAdam/m/3/kernel
s
#Adam/m/3/kernel/Read/ReadVariableOpReadVariableOpAdam/m/3/kernel*
_output_shapes

:@ *
dtype0
t
Adam/v/11/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*
shared_nameAdam/v/11/bias
m
"Adam/v/11/bias/Read/ReadVariableOpReadVariableOpAdam/v/11/bias*
_output_shapes
:@*
dtype0
t
Adam/m/11/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*
shared_nameAdam/m/11/bias
m
"Adam/m/11/bias/Read/ReadVariableOpReadVariableOpAdam/m/11/bias*
_output_shapes
:@*
dtype0
|
Adam/v/11/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@@*!
shared_nameAdam/v/11/kernel
u
$Adam/v/11/kernel/Read/ReadVariableOpReadVariableOpAdam/v/11/kernel*
_output_shapes

:@@*
dtype0
|
Adam/m/11/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@@*!
shared_nameAdam/m/11/kernel
u
$Adam/m/11/kernel/Read/ReadVariableOpReadVariableOpAdam/m/11/kernel*
_output_shapes

:@@*
dtype0
r
Adam/v/8/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*
shared_nameAdam/v/8/bias
k
!Adam/v/8/bias/Read/ReadVariableOpReadVariableOpAdam/v/8/bias*
_output_shapes
:@*
dtype0
r
Adam/m/8/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*
shared_nameAdam/m/8/bias
k
!Adam/m/8/bias/Read/ReadVariableOpReadVariableOpAdam/m/8/bias*
_output_shapes
:@*
dtype0
z
Adam/v/8/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@@* 
shared_nameAdam/v/8/kernel
s
#Adam/v/8/kernel/Read/ReadVariableOpReadVariableOpAdam/v/8/kernel*
_output_shapes

:@@*
dtype0
z
Adam/m/8/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@@* 
shared_nameAdam/m/8/kernel
s
#Adam/m/8/kernel/Read/ReadVariableOpReadVariableOpAdam/m/8/kernel*
_output_shapes

:@@*
dtype0
r
Adam/v/5/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*
shared_nameAdam/v/5/bias
k
!Adam/v/5/bias/Read/ReadVariableOpReadVariableOpAdam/v/5/bias*
_output_shapes
:@*
dtype0
r
Adam/m/5/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*
shared_nameAdam/m/5/bias
k
!Adam/m/5/bias/Read/ReadVariableOpReadVariableOpAdam/m/5/bias*
_output_shapes
:@*
dtype0
z
Adam/v/5/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@@* 
shared_nameAdam/v/5/kernel
s
#Adam/v/5/kernel/Read/ReadVariableOpReadVariableOpAdam/v/5/kernel*
_output_shapes

:@@*
dtype0
z
Adam/m/5/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@@* 
shared_nameAdam/m/5/kernel
s
#Adam/m/5/kernel/Read/ReadVariableOpReadVariableOpAdam/m/5/kernel*
_output_shapes

:@@*
dtype0
r
Adam/v/2/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*
shared_nameAdam/v/2/bias
k
!Adam/v/2/bias/Read/ReadVariableOpReadVariableOpAdam/v/2/bias*
_output_shapes
:@*
dtype0
r
Adam/m/2/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*
shared_nameAdam/m/2/bias
k
!Adam/m/2/bias/Read/ReadVariableOpReadVariableOpAdam/m/2/bias*
_output_shapes
:@*
dtype0
z
Adam/v/2/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@@* 
shared_nameAdam/v/2/kernel
s
#Adam/v/2/kernel/Read/ReadVariableOpReadVariableOpAdam/v/2/kernel*
_output_shapes

:@@*
dtype0
z
Adam/m/2/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@@* 
shared_nameAdam/m/2/kernel
s
#Adam/m/2/kernel/Read/ReadVariableOpReadVariableOpAdam/m/2/kernel*
_output_shapes

:@@*
dtype0
t
Adam/v/10/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*
shared_nameAdam/v/10/bias
m
"Adam/v/10/bias/Read/ReadVariableOpReadVariableOpAdam/v/10/bias*
_output_shapes
:@*
dtype0
t
Adam/m/10/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*
shared_nameAdam/m/10/bias
m
"Adam/m/10/bias/Read/ReadVariableOpReadVariableOpAdam/m/10/bias*
_output_shapes
:@*
dtype0
|
Adam/v/10/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@@*!
shared_nameAdam/v/10/kernel
u
$Adam/v/10/kernel/Read/ReadVariableOpReadVariableOpAdam/v/10/kernel*
_output_shapes

:@@*
dtype0
|
Adam/m/10/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@@*!
shared_nameAdam/m/10/kernel
u
$Adam/m/10/kernel/Read/ReadVariableOpReadVariableOpAdam/m/10/kernel*
_output_shapes

:@@*
dtype0
r
Adam/v/7/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*
shared_nameAdam/v/7/bias
k
!Adam/v/7/bias/Read/ReadVariableOpReadVariableOpAdam/v/7/bias*
_output_shapes
:@*
dtype0
r
Adam/m/7/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*
shared_nameAdam/m/7/bias
k
!Adam/m/7/bias/Read/ReadVariableOpReadVariableOpAdam/m/7/bias*
_output_shapes
:@*
dtype0
z
Adam/v/7/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@@* 
shared_nameAdam/v/7/kernel
s
#Adam/v/7/kernel/Read/ReadVariableOpReadVariableOpAdam/v/7/kernel*
_output_shapes

:@@*
dtype0
z
Adam/m/7/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@@* 
shared_nameAdam/m/7/kernel
s
#Adam/m/7/kernel/Read/ReadVariableOpReadVariableOpAdam/m/7/kernel*
_output_shapes

:@@*
dtype0
r
Adam/v/4/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*
shared_nameAdam/v/4/bias
k
!Adam/v/4/bias/Read/ReadVariableOpReadVariableOpAdam/v/4/bias*
_output_shapes
:@*
dtype0
r
Adam/m/4/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*
shared_nameAdam/m/4/bias
k
!Adam/m/4/bias/Read/ReadVariableOpReadVariableOpAdam/m/4/bias*
_output_shapes
:@*
dtype0
z
Adam/v/4/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@@* 
shared_nameAdam/v/4/kernel
s
#Adam/v/4/kernel/Read/ReadVariableOpReadVariableOpAdam/v/4/kernel*
_output_shapes

:@@*
dtype0
z
Adam/m/4/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@@* 
shared_nameAdam/m/4/kernel
s
#Adam/m/4/kernel/Read/ReadVariableOpReadVariableOpAdam/m/4/kernel*
_output_shapes

:@@*
dtype0
r
Adam/v/1/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*
shared_nameAdam/v/1/bias
k
!Adam/v/1/bias/Read/ReadVariableOpReadVariableOpAdam/v/1/bias*
_output_shapes
:@*
dtype0
r
Adam/m/1/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*
shared_nameAdam/m/1/bias
k
!Adam/m/1/bias/Read/ReadVariableOpReadVariableOpAdam/m/1/bias*
_output_shapes
:@*
dtype0
z
Adam/v/1/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@@* 
shared_nameAdam/v/1/kernel
s
#Adam/v/1/kernel/Read/ReadVariableOpReadVariableOpAdam/v/1/kernel*
_output_shapes

:@@*
dtype0
z
Adam/m/1/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@@* 
shared_nameAdam/m/1/kernel
s
#Adam/m/1/kernel/Read/ReadVariableOpReadVariableOpAdam/m/1/kernel*
_output_shapes

:@@*
dtype0
r
Adam/v/x/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*
shared_nameAdam/v/x/bias
k
!Adam/v/x/bias/Read/ReadVariableOpReadVariableOpAdam/v/x/bias*
_output_shapes
:@*
dtype0
r
Adam/m/x/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*
shared_nameAdam/m/x/bias
k
!Adam/m/x/bias/Read/ReadVariableOpReadVariableOpAdam/m/x/bias*
_output_shapes
:@*
dtype0
z
Adam/v/x/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@@* 
shared_nameAdam/v/x/kernel
s
#Adam/v/x/kernel/Read/ReadVariableOpReadVariableOpAdam/v/x/kernel*
_output_shapes

:@@*
dtype0
z
Adam/m/x/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@@* 
shared_nameAdam/m/x/kernel
s
#Adam/m/x/kernel/Read/ReadVariableOpReadVariableOpAdam/m/x/kernel*
_output_shapes

:@@*
dtype0
t
Adam/v/l3/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*
shared_nameAdam/v/l3/bias
m
"Adam/v/l3/bias/Read/ReadVariableOpReadVariableOpAdam/v/l3/bias*
_output_shapes
:@*
dtype0
t
Adam/m/l3/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*
shared_nameAdam/m/l3/bias
m
"Adam/m/l3/bias/Read/ReadVariableOpReadVariableOpAdam/m/l3/bias*
_output_shapes
:@*
dtype0
|
Adam/v/l3/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@@*!
shared_nameAdam/v/l3/kernel
u
$Adam/v/l3/kernel/Read/ReadVariableOpReadVariableOpAdam/v/l3/kernel*
_output_shapes

:@@*
dtype0
|
Adam/m/l3/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@@*!
shared_nameAdam/m/l3/kernel
u
$Adam/m/l3/kernel/Read/ReadVariableOpReadVariableOpAdam/m/l3/kernel*
_output_shapes

:@@*
dtype0
t
Adam/v/l2/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*
shared_nameAdam/v/l2/bias
m
"Adam/v/l2/bias/Read/ReadVariableOpReadVariableOpAdam/v/l2/bias*
_output_shapes
:@*
dtype0
t
Adam/m/l2/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*
shared_nameAdam/m/l2/bias
m
"Adam/m/l2/bias/Read/ReadVariableOpReadVariableOpAdam/m/l2/bias*
_output_shapes
:@*
dtype0
|
Adam/v/l2/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@*!
shared_nameAdam/v/l2/kernel
u
$Adam/v/l2/kernel/Read/ReadVariableOpReadVariableOpAdam/v/l2/kernel*
_output_shapes

:@*
dtype0
|
Adam/m/l2/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@*!
shared_nameAdam/m/l2/kernel
u
$Adam/m/l2/kernel/Read/ReadVariableOpReadVariableOpAdam/m/l2/kernel*
_output_shapes

:@*
dtype0
n
learning_rateVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_namelearning_rate
g
!learning_rate/Read/ReadVariableOpReadVariableOplearning_rate*
_output_shapes
: *
dtype0
f
	iterationVarHandleOp*
_output_shapes
: *
dtype0	*
shape: *
shared_name	iteration
_
iteration/Read/ReadVariableOpReadVariableOp	iteration*
_output_shapes
: *
dtype0	
d
c/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_namec/bias
]
c/bias/Read/ReadVariableOpReadVariableOpc/bias*
_output_shapes
:*
dtype0
l
c/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
: *
shared_name
c/kernel
e
c/kernel/Read/ReadVariableOpReadVariableOpc/kernel*
_output_shapes

: *
dtype0
h
phi/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_name
phi/bias
a
phi/bias/Read/ReadVariableOpReadVariableOpphi/bias*
_output_shapes
:*
dtype0
p

phi/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
: *
shared_name
phi/kernel
i
phi/kernel/Read/ReadVariableOpReadVariableOp
phi/kernel*
_output_shapes

: *
dtype0
d
P/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_nameP/bias
]
P/bias/Read/ReadVariableOpReadVariableOpP/bias*
_output_shapes
:*
dtype0
l
P/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
: *
shared_name
P/kernel
e
P/kernel/Read/ReadVariableOpReadVariableOpP/kernel*
_output_shapes

: *
dtype0
d
T/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_nameT/bias
]
T/bias/Read/ReadVariableOpReadVariableOpT/bias*
_output_shapes
:*
dtype0
l
T/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
: *
shared_name
T/kernel
e
T/kernel/Read/ReadVariableOpReadVariableOpT/kernel*
_output_shapes

: *
dtype0
f
12/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name	12/bias
_
12/bias/Read/ReadVariableOpReadVariableOp12/bias*
_output_shapes
: *
dtype0
n
	12/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@ *
shared_name	12/kernel
g
12/kernel/Read/ReadVariableOpReadVariableOp	12/kernel*
_output_shapes

:@ *
dtype0
d
9/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name9/bias
]
9/bias/Read/ReadVariableOpReadVariableOp9/bias*
_output_shapes
: *
dtype0
l
9/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@ *
shared_name
9/kernel
e
9/kernel/Read/ReadVariableOpReadVariableOp9/kernel*
_output_shapes

:@ *
dtype0
d
6/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name6/bias
]
6/bias/Read/ReadVariableOpReadVariableOp6/bias*
_output_shapes
: *
dtype0
l
6/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@ *
shared_name
6/kernel
e
6/kernel/Read/ReadVariableOpReadVariableOp6/kernel*
_output_shapes

:@ *
dtype0
d
3/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape: *
shared_name3/bias
]
3/bias/Read/ReadVariableOpReadVariableOp3/bias*
_output_shapes
: *
dtype0
l
3/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@ *
shared_name
3/kernel
e
3/kernel/Read/ReadVariableOpReadVariableOp3/kernel*
_output_shapes

:@ *
dtype0
f
11/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*
shared_name	11/bias
_
11/bias/Read/ReadVariableOpReadVariableOp11/bias*
_output_shapes
:@*
dtype0
n
	11/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@@*
shared_name	11/kernel
g
11/kernel/Read/ReadVariableOpReadVariableOp	11/kernel*
_output_shapes

:@@*
dtype0
d
8/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*
shared_name8/bias
]
8/bias/Read/ReadVariableOpReadVariableOp8/bias*
_output_shapes
:@*
dtype0
l
8/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@@*
shared_name
8/kernel
e
8/kernel/Read/ReadVariableOpReadVariableOp8/kernel*
_output_shapes

:@@*
dtype0
d
5/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*
shared_name5/bias
]
5/bias/Read/ReadVariableOpReadVariableOp5/bias*
_output_shapes
:@*
dtype0
l
5/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@@*
shared_name
5/kernel
e
5/kernel/Read/ReadVariableOpReadVariableOp5/kernel*
_output_shapes

:@@*
dtype0
d
2/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*
shared_name2/bias
]
2/bias/Read/ReadVariableOpReadVariableOp2/bias*
_output_shapes
:@*
dtype0
l
2/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@@*
shared_name
2/kernel
e
2/kernel/Read/ReadVariableOpReadVariableOp2/kernel*
_output_shapes

:@@*
dtype0
f
10/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*
shared_name	10/bias
_
10/bias/Read/ReadVariableOpReadVariableOp10/bias*
_output_shapes
:@*
dtype0
n
	10/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@@*
shared_name	10/kernel
g
10/kernel/Read/ReadVariableOpReadVariableOp	10/kernel*
_output_shapes

:@@*
dtype0
d
7/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*
shared_name7/bias
]
7/bias/Read/ReadVariableOpReadVariableOp7/bias*
_output_shapes
:@*
dtype0
l
7/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@@*
shared_name
7/kernel
e
7/kernel/Read/ReadVariableOpReadVariableOp7/kernel*
_output_shapes

:@@*
dtype0
d
4/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*
shared_name4/bias
]
4/bias/Read/ReadVariableOpReadVariableOp4/bias*
_output_shapes
:@*
dtype0
l
4/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@@*
shared_name
4/kernel
e
4/kernel/Read/ReadVariableOpReadVariableOp4/kernel*
_output_shapes

:@@*
dtype0
d
1/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*
shared_name1/bias
]
1/bias/Read/ReadVariableOpReadVariableOp1/bias*
_output_shapes
:@*
dtype0
l
1/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@@*
shared_name
1/kernel
e
1/kernel/Read/ReadVariableOpReadVariableOp1/kernel*
_output_shapes

:@@*
dtype0
d
x/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*
shared_namex/bias
]
x/bias/Read/ReadVariableOpReadVariableOpx/bias*
_output_shapes
:@*
dtype0
l
x/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@@*
shared_name
x/kernel
e
x/kernel/Read/ReadVariableOpReadVariableOpx/kernel*
_output_shapes

:@@*
dtype0
f
l3/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*
shared_name	l3/bias
_
l3/bias/Read/ReadVariableOpReadVariableOpl3/bias*
_output_shapes
:@*
dtype0
n
	l3/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@@*
shared_name	l3/kernel
g
l3/kernel/Read/ReadVariableOpReadVariableOp	l3/kernel*
_output_shapes

:@@*
dtype0
f
l2/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*
shared_name	l2/bias
_
l2/bias/Read/ReadVariableOpReadVariableOpl2/bias*
_output_shapes
:@*
dtype0
n
	l2/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape
:@*
shared_name	l2/kernel
g
l2/kernel/Read/ReadVariableOpReadVariableOp	l2/kernel*
_output_shapes

:@*
dtype0
b
count_5VarHandleOp*
_output_shapes
: *
dtype0	*
shape: *
shared_name	count_5
[
count_5/Read/ReadVariableOpReadVariableOpcount_5*
_output_shapes
: *
dtype0	
h
varianceVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_name
variance
a
variance/Read/ReadVariableOpReadVariableOpvariance*
_output_shapes
:*
dtype0
`
meanVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_namemean
Y
mean/Read/ReadVariableOpReadVariableOpmean*
_output_shapes
:*
dtype0
x
serving_default_inputPlaceholder*'
_output_shapes
:���������*
dtype0*
shape:���������
�
StatefulPartitionedCallStatefulPartitionedCallserving_default_inputConst_1Const	l2/kernell2/bias	l3/kernell3/biasx/kernelx/bias	10/kernel10/bias7/kernel7/bias4/kernel4/bias1/kernel1/bias	11/kernel11/bias8/kernel8/bias5/kernel5/bias2/kernel2/bias	12/kernel12/bias9/kernel9/bias6/kernel6/bias3/kernel3/biasc/kernelc/bias
phi/kernelphi/biasP/kernelP/biasT/kernelT/bias*4
Tin-
+2)*
Tout
2*
_collective_manager_ids
 *`
_output_shapesN
L:���������:���������:���������:���������*H
_read_only_resource_inputs*
(&	
 !"#$%&'(*-
config_proto

CPU

GPU 2J 8� *.
f)R'
%__inference_signature_wrapper_2218179

NoOpNoOp
��
Const_2Const"/device:CPU:0*
_output_shapes
: *
dtype0*��
value��B�� B��
�
layer-0
layer_with_weights-0
layer-1
layer_with_weights-1
layer-2
layer_with_weights-2
layer-3
layer_with_weights-3
layer-4
layer_with_weights-4
layer-5
layer_with_weights-5
layer-6
layer_with_weights-6
layer-7
	layer_with_weights-7
	layer-8

layer_with_weights-8

layer-9
layer_with_weights-9
layer-10
layer_with_weights-10
layer-11
layer_with_weights-11
layer-12
layer_with_weights-12
layer-13
layer_with_weights-13
layer-14
layer_with_weights-14
layer-15
layer_with_weights-15
layer-16
layer_with_weights-16
layer-17
layer_with_weights-17
layer-18
layer_with_weights-18
layer-19
layer_with_weights-19
layer-20
	variables
trainable_variables
regularization_losses
	keras_api
__call__
*&call_and_return_all_conditional_losses
_default_save_signature
	optimizer
loss

signatures*
* 
�
 	keras_api
!
_keep_axis
"_reduce_axis
#_reduce_axis_mask
$_broadcast_shape
%mean
%
adapt_mean
&variance
&adapt_variance
	'count
(_adapt_function*
�
)	variables
*trainable_variables
+regularization_losses
,	keras_api
-__call__
*.&call_and_return_all_conditional_losses

/kernel
0bias*
�
1	variables
2trainable_variables
3regularization_losses
4	keras_api
5__call__
*6&call_and_return_all_conditional_losses

7kernel
8bias*
�
9	variables
:trainable_variables
;regularization_losses
<	keras_api
=__call__
*>&call_and_return_all_conditional_losses

?kernel
@bias*
�
A	variables
Btrainable_variables
Cregularization_losses
D	keras_api
E__call__
*F&call_and_return_all_conditional_losses

Gkernel
Hbias*
�
I	variables
Jtrainable_variables
Kregularization_losses
L	keras_api
M__call__
*N&call_and_return_all_conditional_losses

Okernel
Pbias*
�
Q	variables
Rtrainable_variables
Sregularization_losses
T	keras_api
U__call__
*V&call_and_return_all_conditional_losses

Wkernel
Xbias*
�
Y	variables
Ztrainable_variables
[regularization_losses
\	keras_api
]__call__
*^&call_and_return_all_conditional_losses

_kernel
`bias*
�
a	variables
btrainable_variables
cregularization_losses
d	keras_api
e__call__
*f&call_and_return_all_conditional_losses

gkernel
hbias*
�
i	variables
jtrainable_variables
kregularization_losses
l	keras_api
m__call__
*n&call_and_return_all_conditional_losses

okernel
pbias*
�
q	variables
rtrainable_variables
sregularization_losses
t	keras_api
u__call__
*v&call_and_return_all_conditional_losses

wkernel
xbias*
�
y	variables
ztrainable_variables
{regularization_losses
|	keras_api
}__call__
*~&call_and_return_all_conditional_losses

kernel
	�bias*
�
�	variables
�trainable_variables
�regularization_losses
�	keras_api
�__call__
+�&call_and_return_all_conditional_losses
�kernel
	�bias*
�
�	variables
�trainable_variables
�regularization_losses
�	keras_api
�__call__
+�&call_and_return_all_conditional_losses
�kernel
	�bias*
�
�	variables
�trainable_variables
�regularization_losses
�	keras_api
�__call__
+�&call_and_return_all_conditional_losses
�kernel
	�bias*
�
�	variables
�trainable_variables
�regularization_losses
�	keras_api
�__call__
+�&call_and_return_all_conditional_losses
�kernel
	�bias*
�
�	variables
�trainable_variables
�regularization_losses
�	keras_api
�__call__
+�&call_and_return_all_conditional_losses
�kernel
	�bias*
�
�	variables
�trainable_variables
�regularization_losses
�	keras_api
�__call__
+�&call_and_return_all_conditional_losses
�kernel
	�bias*
�
�	variables
�trainable_variables
�regularization_losses
�	keras_api
�__call__
+�&call_and_return_all_conditional_losses
�kernel
	�bias*
�
�	variables
�trainable_variables
�regularization_losses
�	keras_api
�__call__
+�&call_and_return_all_conditional_losses
�kernel
	�bias*
�
%0
&1
'2
/3
04
75
86
?7
@8
G9
H10
O11
P12
W13
X14
_15
`16
g17
h18
o19
p20
w21
x22
23
�24
�25
�26
�27
�28
�29
�30
�31
�32
�33
�34
�35
�36
�37
�38
�39
�40*
�
/0
01
72
83
?4
@5
G6
H7
O8
P9
W10
X11
_12
`13
g14
h15
o16
p17
w18
x19
20
�21
�22
�23
�24
�25
�26
�27
�28
�29
�30
�31
�32
�33
�34
�35
�36
�37*
* 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
	variables
trainable_variables
regularization_losses
__call__
_default_save_signature
*&call_and_return_all_conditional_losses
&"call_and_return_conditional_losses*
:
�trace_0
�trace_1
�trace_2
�trace_3* 
:
�trace_0
�trace_1
�trace_2
�trace_3* 
"
�	capture_0
�	capture_1* 
�
�
_variables
�_iterations
�_learning_rate
�_index_dict
�
_momentums
�_velocities
�_update_step_xla*
* 

�serving_default* 
* 
* 
* 
* 
* 
RL
VARIABLE_VALUEmean4layer_with_weights-0/mean/.ATTRIBUTES/VARIABLE_VALUE*
ZT
VARIABLE_VALUEvariance8layer_with_weights-0/variance/.ATTRIBUTES/VARIABLE_VALUE*
VP
VARIABLE_VALUEcount_55layer_with_weights-0/count/.ATTRIBUTES/VARIABLE_VALUE*

�trace_0* 

/0
01*

/0
01*
* 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
)	variables
*trainable_variables
+regularization_losses
-__call__
*.&call_and_return_all_conditional_losses
&."call_and_return_conditional_losses*

�trace_0* 

�trace_0* 
YS
VARIABLE_VALUE	l2/kernel6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUE*
UO
VARIABLE_VALUEl2/bias4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUE*

70
81*

70
81*
* 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
1	variables
2trainable_variables
3regularization_losses
5__call__
*6&call_and_return_all_conditional_losses
&6"call_and_return_conditional_losses*

�trace_0* 

�trace_0* 
YS
VARIABLE_VALUE	l3/kernel6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUE*
UO
VARIABLE_VALUEl3/bias4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUE*

?0
@1*

?0
@1*
* 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
9	variables
:trainable_variables
;regularization_losses
=__call__
*>&call_and_return_all_conditional_losses
&>"call_and_return_conditional_losses*

�trace_0* 

�trace_0* 
XR
VARIABLE_VALUEx/kernel6layer_with_weights-3/kernel/.ATTRIBUTES/VARIABLE_VALUE*
TN
VARIABLE_VALUEx/bias4layer_with_weights-3/bias/.ATTRIBUTES/VARIABLE_VALUE*

G0
H1*

G0
H1*
* 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
A	variables
Btrainable_variables
Cregularization_losses
E__call__
*F&call_and_return_all_conditional_losses
&F"call_and_return_conditional_losses*

�trace_0* 

�trace_0* 
XR
VARIABLE_VALUE1/kernel6layer_with_weights-4/kernel/.ATTRIBUTES/VARIABLE_VALUE*
TN
VARIABLE_VALUE1/bias4layer_with_weights-4/bias/.ATTRIBUTES/VARIABLE_VALUE*

O0
P1*

O0
P1*
* 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
I	variables
Jtrainable_variables
Kregularization_losses
M__call__
*N&call_and_return_all_conditional_losses
&N"call_and_return_conditional_losses*

�trace_0* 

�trace_0* 
XR
VARIABLE_VALUE4/kernel6layer_with_weights-5/kernel/.ATTRIBUTES/VARIABLE_VALUE*
TN
VARIABLE_VALUE4/bias4layer_with_weights-5/bias/.ATTRIBUTES/VARIABLE_VALUE*

W0
X1*

W0
X1*
* 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
Q	variables
Rtrainable_variables
Sregularization_losses
U__call__
*V&call_and_return_all_conditional_losses
&V"call_and_return_conditional_losses*

�trace_0* 

�trace_0* 
XR
VARIABLE_VALUE7/kernel6layer_with_weights-6/kernel/.ATTRIBUTES/VARIABLE_VALUE*
TN
VARIABLE_VALUE7/bias4layer_with_weights-6/bias/.ATTRIBUTES/VARIABLE_VALUE*

_0
`1*

_0
`1*
* 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
Y	variables
Ztrainable_variables
[regularization_losses
]__call__
*^&call_and_return_all_conditional_losses
&^"call_and_return_conditional_losses*

�trace_0* 

�trace_0* 
YS
VARIABLE_VALUE	10/kernel6layer_with_weights-7/kernel/.ATTRIBUTES/VARIABLE_VALUE*
UO
VARIABLE_VALUE10/bias4layer_with_weights-7/bias/.ATTRIBUTES/VARIABLE_VALUE*

g0
h1*

g0
h1*
* 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
a	variables
btrainable_variables
cregularization_losses
e__call__
*f&call_and_return_all_conditional_losses
&f"call_and_return_conditional_losses*

�trace_0* 

�trace_0* 
XR
VARIABLE_VALUE2/kernel6layer_with_weights-8/kernel/.ATTRIBUTES/VARIABLE_VALUE*
TN
VARIABLE_VALUE2/bias4layer_with_weights-8/bias/.ATTRIBUTES/VARIABLE_VALUE*

o0
p1*

o0
p1*
* 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
i	variables
jtrainable_variables
kregularization_losses
m__call__
*n&call_and_return_all_conditional_losses
&n"call_and_return_conditional_losses*

�trace_0* 

�trace_0* 
XR
VARIABLE_VALUE5/kernel6layer_with_weights-9/kernel/.ATTRIBUTES/VARIABLE_VALUE*
TN
VARIABLE_VALUE5/bias4layer_with_weights-9/bias/.ATTRIBUTES/VARIABLE_VALUE*

w0
x1*

w0
x1*
* 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
q	variables
rtrainable_variables
sregularization_losses
u__call__
*v&call_and_return_all_conditional_losses
&v"call_and_return_conditional_losses*

�trace_0* 

�trace_0* 
YS
VARIABLE_VALUE8/kernel7layer_with_weights-10/kernel/.ATTRIBUTES/VARIABLE_VALUE*
UO
VARIABLE_VALUE8/bias5layer_with_weights-10/bias/.ATTRIBUTES/VARIABLE_VALUE*

0
�1*

0
�1*
* 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
y	variables
ztrainable_variables
{regularization_losses
}__call__
*~&call_and_return_all_conditional_losses
&~"call_and_return_conditional_losses*

�trace_0* 

�trace_0* 
ZT
VARIABLE_VALUE	11/kernel7layer_with_weights-11/kernel/.ATTRIBUTES/VARIABLE_VALUE*
VP
VARIABLE_VALUE11/bias5layer_with_weights-11/bias/.ATTRIBUTES/VARIABLE_VALUE*

�0
�1*

�0
�1*
* 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
�	variables
�trainable_variables
�regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses*

�trace_0* 

�trace_0* 
YS
VARIABLE_VALUE3/kernel7layer_with_weights-12/kernel/.ATTRIBUTES/VARIABLE_VALUE*
UO
VARIABLE_VALUE3/bias5layer_with_weights-12/bias/.ATTRIBUTES/VARIABLE_VALUE*

�0
�1*

�0
�1*
* 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
�	variables
�trainable_variables
�regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses*

�trace_0* 

�trace_0* 
YS
VARIABLE_VALUE6/kernel7layer_with_weights-13/kernel/.ATTRIBUTES/VARIABLE_VALUE*
UO
VARIABLE_VALUE6/bias5layer_with_weights-13/bias/.ATTRIBUTES/VARIABLE_VALUE*

�0
�1*

�0
�1*
* 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
�	variables
�trainable_variables
�regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses*

�trace_0* 

�trace_0* 
YS
VARIABLE_VALUE9/kernel7layer_with_weights-14/kernel/.ATTRIBUTES/VARIABLE_VALUE*
UO
VARIABLE_VALUE9/bias5layer_with_weights-14/bias/.ATTRIBUTES/VARIABLE_VALUE*

�0
�1*

�0
�1*
* 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
�	variables
�trainable_variables
�regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses*

�trace_0* 

�trace_0* 
ZT
VARIABLE_VALUE	12/kernel7layer_with_weights-15/kernel/.ATTRIBUTES/VARIABLE_VALUE*
VP
VARIABLE_VALUE12/bias5layer_with_weights-15/bias/.ATTRIBUTES/VARIABLE_VALUE*

�0
�1*

�0
�1*
* 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
�	variables
�trainable_variables
�regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses*

�trace_0* 

�trace_0* 
YS
VARIABLE_VALUET/kernel7layer_with_weights-16/kernel/.ATTRIBUTES/VARIABLE_VALUE*
UO
VARIABLE_VALUET/bias5layer_with_weights-16/bias/.ATTRIBUTES/VARIABLE_VALUE*

�0
�1*

�0
�1*
* 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
�	variables
�trainable_variables
�regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses*

�trace_0* 

�trace_0* 
YS
VARIABLE_VALUEP/kernel7layer_with_weights-17/kernel/.ATTRIBUTES/VARIABLE_VALUE*
UO
VARIABLE_VALUEP/bias5layer_with_weights-17/bias/.ATTRIBUTES/VARIABLE_VALUE*

�0
�1*

�0
�1*
* 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
�	variables
�trainable_variables
�regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses*

�trace_0* 

�trace_0* 
[U
VARIABLE_VALUE
phi/kernel7layer_with_weights-18/kernel/.ATTRIBUTES/VARIABLE_VALUE*
WQ
VARIABLE_VALUEphi/bias5layer_with_weights-18/bias/.ATTRIBUTES/VARIABLE_VALUE*

�0
�1*

�0
�1*
* 
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
�	variables
�trainable_variables
�regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses*

�trace_0* 

�trace_0* 
YS
VARIABLE_VALUEc/kernel7layer_with_weights-19/kernel/.ATTRIBUTES/VARIABLE_VALUE*
UO
VARIABLE_VALUEc/bias5layer_with_weights-19/bias/.ATTRIBUTES/VARIABLE_VALUE*

%0
&1
'2*
�
0
1
2
3
4
5
6
7
	8

9
10
11
12
13
14
15
16
17
18
19
20*
,
�0
�1
�2
�3
�4*
* 
* 
"
�	capture_0
�	capture_1* 
"
�	capture_0
�	capture_1* 
"
�	capture_0
�	capture_1* 
"
�	capture_0
�	capture_1* 
"
�	capture_0
�	capture_1* 
"
�	capture_0
�	capture_1* 
"
�	capture_0
�	capture_1* 
"
�	capture_0
�	capture_1* 
* 
* 
�
�0
�1
�2
�3
�4
�5
�6
�7
�8
�9
�10
�11
�12
�13
�14
�15
�16
�17
�18
�19
�20
�21
�22
�23
�24
�25
�26
�27
�28
�29
�30
�31
�32
�33
�34
�35
�36
�37
�38
�39
�40
�41
�42
�43
�44
�45
�46
�47
�48
�49
�50
�51
�52
�53
�54
�55
�56
�57
�58
�59
�60
�61
�62
�63
�64
�65
�66
�67
�68
�69
�70
�71
�72
�73
�74
�75
�76*
SM
VARIABLE_VALUE	iteration0optimizer/_iterations/.ATTRIBUTES/VARIABLE_VALUE*
ZT
VARIABLE_VALUElearning_rate3optimizer/_learning_rate/.ATTRIBUTES/VARIABLE_VALUE*
* 
�
�0
�1
�2
�3
�4
�5
�6
�7
�8
�9
�10
�11
�12
�13
�14
�15
�16
�17
�18
�19
�20
�21
�22
�23
�24
�25
�26
�27
�28
�29
�30
�31
�32
�33
�34
�35
�36
�37*
�
�0
�1
�2
�3
�4
�5
�6
�7
�8
�9
�10
�11
�12
�13
�14
�15
�16
�17
�18
�19
�20
�21
�22
�23
�24
�25
�26
�27
�28
�29
�30
�31
�32
�33
�34
�35
�36
�37*
* 
"
�	capture_0
�	capture_1* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
* 
<
�	variables
�	keras_api

�total

�count*
<
�	variables
�	keras_api

�total

�count*
<
�	variables
�	keras_api

�total

�count*
<
�	variables
�	keras_api

�total

�count*
<
�	variables
�	keras_api

�total

�count*
[U
VARIABLE_VALUEAdam/m/l2/kernel1optimizer/_variables/1/.ATTRIBUTES/VARIABLE_VALUE*
[U
VARIABLE_VALUEAdam/v/l2/kernel1optimizer/_variables/2/.ATTRIBUTES/VARIABLE_VALUE*
YS
VARIABLE_VALUEAdam/m/l2/bias1optimizer/_variables/3/.ATTRIBUTES/VARIABLE_VALUE*
YS
VARIABLE_VALUEAdam/v/l2/bias1optimizer/_variables/4/.ATTRIBUTES/VARIABLE_VALUE*
[U
VARIABLE_VALUEAdam/m/l3/kernel1optimizer/_variables/5/.ATTRIBUTES/VARIABLE_VALUE*
[U
VARIABLE_VALUEAdam/v/l3/kernel1optimizer/_variables/6/.ATTRIBUTES/VARIABLE_VALUE*
YS
VARIABLE_VALUEAdam/m/l3/bias1optimizer/_variables/7/.ATTRIBUTES/VARIABLE_VALUE*
YS
VARIABLE_VALUEAdam/v/l3/bias1optimizer/_variables/8/.ATTRIBUTES/VARIABLE_VALUE*
ZT
VARIABLE_VALUEAdam/m/x/kernel1optimizer/_variables/9/.ATTRIBUTES/VARIABLE_VALUE*
[U
VARIABLE_VALUEAdam/v/x/kernel2optimizer/_variables/10/.ATTRIBUTES/VARIABLE_VALUE*
YS
VARIABLE_VALUEAdam/m/x/bias2optimizer/_variables/11/.ATTRIBUTES/VARIABLE_VALUE*
YS
VARIABLE_VALUEAdam/v/x/bias2optimizer/_variables/12/.ATTRIBUTES/VARIABLE_VALUE*
[U
VARIABLE_VALUEAdam/m/1/kernel2optimizer/_variables/13/.ATTRIBUTES/VARIABLE_VALUE*
[U
VARIABLE_VALUEAdam/v/1/kernel2optimizer/_variables/14/.ATTRIBUTES/VARIABLE_VALUE*
YS
VARIABLE_VALUEAdam/m/1/bias2optimizer/_variables/15/.ATTRIBUTES/VARIABLE_VALUE*
YS
VARIABLE_VALUEAdam/v/1/bias2optimizer/_variables/16/.ATTRIBUTES/VARIABLE_VALUE*
[U
VARIABLE_VALUEAdam/m/4/kernel2optimizer/_variables/17/.ATTRIBUTES/VARIABLE_VALUE*
[U
VARIABLE_VALUEAdam/v/4/kernel2optimizer/_variables/18/.ATTRIBUTES/VARIABLE_VALUE*
YS
VARIABLE_VALUEAdam/m/4/bias2optimizer/_variables/19/.ATTRIBUTES/VARIABLE_VALUE*
YS
VARIABLE_VALUEAdam/v/4/bias2optimizer/_variables/20/.ATTRIBUTES/VARIABLE_VALUE*
[U
VARIABLE_VALUEAdam/m/7/kernel2optimizer/_variables/21/.ATTRIBUTES/VARIABLE_VALUE*
[U
VARIABLE_VALUEAdam/v/7/kernel2optimizer/_variables/22/.ATTRIBUTES/VARIABLE_VALUE*
YS
VARIABLE_VALUEAdam/m/7/bias2optimizer/_variables/23/.ATTRIBUTES/VARIABLE_VALUE*
YS
VARIABLE_VALUEAdam/v/7/bias2optimizer/_variables/24/.ATTRIBUTES/VARIABLE_VALUE*
\V
VARIABLE_VALUEAdam/m/10/kernel2optimizer/_variables/25/.ATTRIBUTES/VARIABLE_VALUE*
\V
VARIABLE_VALUEAdam/v/10/kernel2optimizer/_variables/26/.ATTRIBUTES/VARIABLE_VALUE*
ZT
VARIABLE_VALUEAdam/m/10/bias2optimizer/_variables/27/.ATTRIBUTES/VARIABLE_VALUE*
ZT
VARIABLE_VALUEAdam/v/10/bias2optimizer/_variables/28/.ATTRIBUTES/VARIABLE_VALUE*
[U
VARIABLE_VALUEAdam/m/2/kernel2optimizer/_variables/29/.ATTRIBUTES/VARIABLE_VALUE*
[U
VARIABLE_VALUEAdam/v/2/kernel2optimizer/_variables/30/.ATTRIBUTES/VARIABLE_VALUE*
YS
VARIABLE_VALUEAdam/m/2/bias2optimizer/_variables/31/.ATTRIBUTES/VARIABLE_VALUE*
YS
VARIABLE_VALUEAdam/v/2/bias2optimizer/_variables/32/.ATTRIBUTES/VARIABLE_VALUE*
[U
VARIABLE_VALUEAdam/m/5/kernel2optimizer/_variables/33/.ATTRIBUTES/VARIABLE_VALUE*
[U
VARIABLE_VALUEAdam/v/5/kernel2optimizer/_variables/34/.ATTRIBUTES/VARIABLE_VALUE*
YS
VARIABLE_VALUEAdam/m/5/bias2optimizer/_variables/35/.ATTRIBUTES/VARIABLE_VALUE*
YS
VARIABLE_VALUEAdam/v/5/bias2optimizer/_variables/36/.ATTRIBUTES/VARIABLE_VALUE*
[U
VARIABLE_VALUEAdam/m/8/kernel2optimizer/_variables/37/.ATTRIBUTES/VARIABLE_VALUE*
[U
VARIABLE_VALUEAdam/v/8/kernel2optimizer/_variables/38/.ATTRIBUTES/VARIABLE_VALUE*
YS
VARIABLE_VALUEAdam/m/8/bias2optimizer/_variables/39/.ATTRIBUTES/VARIABLE_VALUE*
YS
VARIABLE_VALUEAdam/v/8/bias2optimizer/_variables/40/.ATTRIBUTES/VARIABLE_VALUE*
\V
VARIABLE_VALUEAdam/m/11/kernel2optimizer/_variables/41/.ATTRIBUTES/VARIABLE_VALUE*
\V
VARIABLE_VALUEAdam/v/11/kernel2optimizer/_variables/42/.ATTRIBUTES/VARIABLE_VALUE*
ZT
VARIABLE_VALUEAdam/m/11/bias2optimizer/_variables/43/.ATTRIBUTES/VARIABLE_VALUE*
ZT
VARIABLE_VALUEAdam/v/11/bias2optimizer/_variables/44/.ATTRIBUTES/VARIABLE_VALUE*
[U
VARIABLE_VALUEAdam/m/3/kernel2optimizer/_variables/45/.ATTRIBUTES/VARIABLE_VALUE*
[U
VARIABLE_VALUEAdam/v/3/kernel2optimizer/_variables/46/.ATTRIBUTES/VARIABLE_VALUE*
YS
VARIABLE_VALUEAdam/m/3/bias2optimizer/_variables/47/.ATTRIBUTES/VARIABLE_VALUE*
YS
VARIABLE_VALUEAdam/v/3/bias2optimizer/_variables/48/.ATTRIBUTES/VARIABLE_VALUE*
[U
VARIABLE_VALUEAdam/m/6/kernel2optimizer/_variables/49/.ATTRIBUTES/VARIABLE_VALUE*
[U
VARIABLE_VALUEAdam/v/6/kernel2optimizer/_variables/50/.ATTRIBUTES/VARIABLE_VALUE*
YS
VARIABLE_VALUEAdam/m/6/bias2optimizer/_variables/51/.ATTRIBUTES/VARIABLE_VALUE*
YS
VARIABLE_VALUEAdam/v/6/bias2optimizer/_variables/52/.ATTRIBUTES/VARIABLE_VALUE*
[U
VARIABLE_VALUEAdam/m/9/kernel2optimizer/_variables/53/.ATTRIBUTES/VARIABLE_VALUE*
[U
VARIABLE_VALUEAdam/v/9/kernel2optimizer/_variables/54/.ATTRIBUTES/VARIABLE_VALUE*
YS
VARIABLE_VALUEAdam/m/9/bias2optimizer/_variables/55/.ATTRIBUTES/VARIABLE_VALUE*
YS
VARIABLE_VALUEAdam/v/9/bias2optimizer/_variables/56/.ATTRIBUTES/VARIABLE_VALUE*
\V
VARIABLE_VALUEAdam/m/12/kernel2optimizer/_variables/57/.ATTRIBUTES/VARIABLE_VALUE*
\V
VARIABLE_VALUEAdam/v/12/kernel2optimizer/_variables/58/.ATTRIBUTES/VARIABLE_VALUE*
ZT
VARIABLE_VALUEAdam/m/12/bias2optimizer/_variables/59/.ATTRIBUTES/VARIABLE_VALUE*
ZT
VARIABLE_VALUEAdam/v/12/bias2optimizer/_variables/60/.ATTRIBUTES/VARIABLE_VALUE*
[U
VARIABLE_VALUEAdam/m/T/kernel2optimizer/_variables/61/.ATTRIBUTES/VARIABLE_VALUE*
[U
VARIABLE_VALUEAdam/v/T/kernel2optimizer/_variables/62/.ATTRIBUTES/VARIABLE_VALUE*
YS
VARIABLE_VALUEAdam/m/T/bias2optimizer/_variables/63/.ATTRIBUTES/VARIABLE_VALUE*
YS
VARIABLE_VALUEAdam/v/T/bias2optimizer/_variables/64/.ATTRIBUTES/VARIABLE_VALUE*
[U
VARIABLE_VALUEAdam/m/P/kernel2optimizer/_variables/65/.ATTRIBUTES/VARIABLE_VALUE*
[U
VARIABLE_VALUEAdam/v/P/kernel2optimizer/_variables/66/.ATTRIBUTES/VARIABLE_VALUE*
YS
VARIABLE_VALUEAdam/m/P/bias2optimizer/_variables/67/.ATTRIBUTES/VARIABLE_VALUE*
YS
VARIABLE_VALUEAdam/v/P/bias2optimizer/_variables/68/.ATTRIBUTES/VARIABLE_VALUE*
]W
VARIABLE_VALUEAdam/m/phi/kernel2optimizer/_variables/69/.ATTRIBUTES/VARIABLE_VALUE*
]W
VARIABLE_VALUEAdam/v/phi/kernel2optimizer/_variables/70/.ATTRIBUTES/VARIABLE_VALUE*
[U
VARIABLE_VALUEAdam/m/phi/bias2optimizer/_variables/71/.ATTRIBUTES/VARIABLE_VALUE*
[U
VARIABLE_VALUEAdam/v/phi/bias2optimizer/_variables/72/.ATTRIBUTES/VARIABLE_VALUE*
[U
VARIABLE_VALUEAdam/m/c/kernel2optimizer/_variables/73/.ATTRIBUTES/VARIABLE_VALUE*
[U
VARIABLE_VALUEAdam/v/c/kernel2optimizer/_variables/74/.ATTRIBUTES/VARIABLE_VALUE*
YS
VARIABLE_VALUEAdam/m/c/bias2optimizer/_variables/75/.ATTRIBUTES/VARIABLE_VALUE*
YS
VARIABLE_VALUEAdam/v/c/bias2optimizer/_variables/76/.ATTRIBUTES/VARIABLE_VALUE*

�0
�1*

�	variables*
UO
VARIABLE_VALUEtotal_44keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUE*
UO
VARIABLE_VALUEcount_44keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUE*

�0
�1*

�	variables*
UO
VARIABLE_VALUEtotal_34keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUE*
UO
VARIABLE_VALUEcount_34keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUE*

�0
�1*

�	variables*
UO
VARIABLE_VALUEtotal_24keras_api/metrics/2/total/.ATTRIBUTES/VARIABLE_VALUE*
UO
VARIABLE_VALUEcount_24keras_api/metrics/2/count/.ATTRIBUTES/VARIABLE_VALUE*

�0
�1*

�	variables*
UO
VARIABLE_VALUEtotal_14keras_api/metrics/3/total/.ATTRIBUTES/VARIABLE_VALUE*
UO
VARIABLE_VALUEcount_14keras_api/metrics/3/count/.ATTRIBUTES/VARIABLE_VALUE*

�0
�1*

�	variables*
SM
VARIABLE_VALUEtotal4keras_api/metrics/4/total/.ATTRIBUTES/VARIABLE_VALUE*
SM
VARIABLE_VALUEcount4keras_api/metrics/4/count/.ATTRIBUTES/VARIABLE_VALUE*
O
saver_filenamePlaceholder*
_output_shapes
: *
dtype0*
shape: 
�%
StatefulPartitionedCall_1StatefulPartitionedCallsaver_filenamemean/Read/ReadVariableOpvariance/Read/ReadVariableOpcount_5/Read/ReadVariableOpl2/kernel/Read/ReadVariableOpl2/bias/Read/ReadVariableOpl3/kernel/Read/ReadVariableOpl3/bias/Read/ReadVariableOpx/kernel/Read/ReadVariableOpx/bias/Read/ReadVariableOp1/kernel/Read/ReadVariableOp1/bias/Read/ReadVariableOp4/kernel/Read/ReadVariableOp4/bias/Read/ReadVariableOp7/kernel/Read/ReadVariableOp7/bias/Read/ReadVariableOp10/kernel/Read/ReadVariableOp10/bias/Read/ReadVariableOp2/kernel/Read/ReadVariableOp2/bias/Read/ReadVariableOp5/kernel/Read/ReadVariableOp5/bias/Read/ReadVariableOp8/kernel/Read/ReadVariableOp8/bias/Read/ReadVariableOp11/kernel/Read/ReadVariableOp11/bias/Read/ReadVariableOp3/kernel/Read/ReadVariableOp3/bias/Read/ReadVariableOp6/kernel/Read/ReadVariableOp6/bias/Read/ReadVariableOp9/kernel/Read/ReadVariableOp9/bias/Read/ReadVariableOp12/kernel/Read/ReadVariableOp12/bias/Read/ReadVariableOpT/kernel/Read/ReadVariableOpT/bias/Read/ReadVariableOpP/kernel/Read/ReadVariableOpP/bias/Read/ReadVariableOpphi/kernel/Read/ReadVariableOpphi/bias/Read/ReadVariableOpc/kernel/Read/ReadVariableOpc/bias/Read/ReadVariableOpiteration/Read/ReadVariableOp!learning_rate/Read/ReadVariableOp$Adam/m/l2/kernel/Read/ReadVariableOp$Adam/v/l2/kernel/Read/ReadVariableOp"Adam/m/l2/bias/Read/ReadVariableOp"Adam/v/l2/bias/Read/ReadVariableOp$Adam/m/l3/kernel/Read/ReadVariableOp$Adam/v/l3/kernel/Read/ReadVariableOp"Adam/m/l3/bias/Read/ReadVariableOp"Adam/v/l3/bias/Read/ReadVariableOp#Adam/m/x/kernel/Read/ReadVariableOp#Adam/v/x/kernel/Read/ReadVariableOp!Adam/m/x/bias/Read/ReadVariableOp!Adam/v/x/bias/Read/ReadVariableOp#Adam/m/1/kernel/Read/ReadVariableOp#Adam/v/1/kernel/Read/ReadVariableOp!Adam/m/1/bias/Read/ReadVariableOp!Adam/v/1/bias/Read/ReadVariableOp#Adam/m/4/kernel/Read/ReadVariableOp#Adam/v/4/kernel/Read/ReadVariableOp!Adam/m/4/bias/Read/ReadVariableOp!Adam/v/4/bias/Read/ReadVariableOp#Adam/m/7/kernel/Read/ReadVariableOp#Adam/v/7/kernel/Read/ReadVariableOp!Adam/m/7/bias/Read/ReadVariableOp!Adam/v/7/bias/Read/ReadVariableOp$Adam/m/10/kernel/Read/ReadVariableOp$Adam/v/10/kernel/Read/ReadVariableOp"Adam/m/10/bias/Read/ReadVariableOp"Adam/v/10/bias/Read/ReadVariableOp#Adam/m/2/kernel/Read/ReadVariableOp#Adam/v/2/kernel/Read/ReadVariableOp!Adam/m/2/bias/Read/ReadVariableOp!Adam/v/2/bias/Read/ReadVariableOp#Adam/m/5/kernel/Read/ReadVariableOp#Adam/v/5/kernel/Read/ReadVariableOp!Adam/m/5/bias/Read/ReadVariableOp!Adam/v/5/bias/Read/ReadVariableOp#Adam/m/8/kernel/Read/ReadVariableOp#Adam/v/8/kernel/Read/ReadVariableOp!Adam/m/8/bias/Read/ReadVariableOp!Adam/v/8/bias/Read/ReadVariableOp$Adam/m/11/kernel/Read/ReadVariableOp$Adam/v/11/kernel/Read/ReadVariableOp"Adam/m/11/bias/Read/ReadVariableOp"Adam/v/11/bias/Read/ReadVariableOp#Adam/m/3/kernel/Read/ReadVariableOp#Adam/v/3/kernel/Read/ReadVariableOp!Adam/m/3/bias/Read/ReadVariableOp!Adam/v/3/bias/Read/ReadVariableOp#Adam/m/6/kernel/Read/ReadVariableOp#Adam/v/6/kernel/Read/ReadVariableOp!Adam/m/6/bias/Read/ReadVariableOp!Adam/v/6/bias/Read/ReadVariableOp#Adam/m/9/kernel/Read/ReadVariableOp#Adam/v/9/kernel/Read/ReadVariableOp!Adam/m/9/bias/Read/ReadVariableOp!Adam/v/9/bias/Read/ReadVariableOp$Adam/m/12/kernel/Read/ReadVariableOp$Adam/v/12/kernel/Read/ReadVariableOp"Adam/m/12/bias/Read/ReadVariableOp"Adam/v/12/bias/Read/ReadVariableOp#Adam/m/T/kernel/Read/ReadVariableOp#Adam/v/T/kernel/Read/ReadVariableOp!Adam/m/T/bias/Read/ReadVariableOp!Adam/v/T/bias/Read/ReadVariableOp#Adam/m/P/kernel/Read/ReadVariableOp#Adam/v/P/kernel/Read/ReadVariableOp!Adam/m/P/bias/Read/ReadVariableOp!Adam/v/P/bias/Read/ReadVariableOp%Adam/m/phi/kernel/Read/ReadVariableOp%Adam/v/phi/kernel/Read/ReadVariableOp#Adam/m/phi/bias/Read/ReadVariableOp#Adam/v/phi/bias/Read/ReadVariableOp#Adam/m/c/kernel/Read/ReadVariableOp#Adam/v/c/kernel/Read/ReadVariableOp!Adam/m/c/bias/Read/ReadVariableOp!Adam/v/c/bias/Read/ReadVariableOptotal_4/Read/ReadVariableOpcount_4/Read/ReadVariableOptotal_3/Read/ReadVariableOpcount_3/Read/ReadVariableOptotal_2/Read/ReadVariableOpcount_2/Read/ReadVariableOptotal_1/Read/ReadVariableOpcount_1/Read/ReadVariableOptotal/Read/ReadVariableOpcount/Read/ReadVariableOpConst_2*�
Tin�
�2�		*
Tout
2*
_collective_manager_ids
 *
_output_shapes
: * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *)
f$R"
 __inference__traced_save_2219495
�
StatefulPartitionedCall_2StatefulPartitionedCallsaver_filenamemeanvariancecount_5	l2/kernell2/bias	l3/kernell3/biasx/kernelx/bias1/kernel1/bias4/kernel4/bias7/kernel7/bias	10/kernel10/bias2/kernel2/bias5/kernel5/bias8/kernel8/bias	11/kernel11/bias3/kernel3/bias6/kernel6/bias9/kernel9/bias	12/kernel12/biasT/kernelT/biasP/kernelP/bias
phi/kernelphi/biasc/kernelc/bias	iterationlearning_rateAdam/m/l2/kernelAdam/v/l2/kernelAdam/m/l2/biasAdam/v/l2/biasAdam/m/l3/kernelAdam/v/l3/kernelAdam/m/l3/biasAdam/v/l3/biasAdam/m/x/kernelAdam/v/x/kernelAdam/m/x/biasAdam/v/x/biasAdam/m/1/kernelAdam/v/1/kernelAdam/m/1/biasAdam/v/1/biasAdam/m/4/kernelAdam/v/4/kernelAdam/m/4/biasAdam/v/4/biasAdam/m/7/kernelAdam/v/7/kernelAdam/m/7/biasAdam/v/7/biasAdam/m/10/kernelAdam/v/10/kernelAdam/m/10/biasAdam/v/10/biasAdam/m/2/kernelAdam/v/2/kernelAdam/m/2/biasAdam/v/2/biasAdam/m/5/kernelAdam/v/5/kernelAdam/m/5/biasAdam/v/5/biasAdam/m/8/kernelAdam/v/8/kernelAdam/m/8/biasAdam/v/8/biasAdam/m/11/kernelAdam/v/11/kernelAdam/m/11/biasAdam/v/11/biasAdam/m/3/kernelAdam/v/3/kernelAdam/m/3/biasAdam/v/3/biasAdam/m/6/kernelAdam/v/6/kernelAdam/m/6/biasAdam/v/6/biasAdam/m/9/kernelAdam/v/9/kernelAdam/m/9/biasAdam/v/9/biasAdam/m/12/kernelAdam/v/12/kernelAdam/m/12/biasAdam/v/12/biasAdam/m/T/kernelAdam/v/T/kernelAdam/m/T/biasAdam/v/T/biasAdam/m/P/kernelAdam/v/P/kernelAdam/m/P/biasAdam/v/P/biasAdam/m/phi/kernelAdam/v/phi/kernelAdam/m/phi/biasAdam/v/phi/biasAdam/m/c/kernelAdam/v/c/kernelAdam/m/c/biasAdam/v/c/biastotal_4count_4total_3count_3total_2count_2total_1count_1totalcount*�
Tin�
�2�*
Tout
2*
_collective_manager_ids
 *
_output_shapes
: * 
_read_only_resource_inputs
 *-
config_proto

CPU

GPU 2J 8� *,
f'R%
#__inference__traced_restore_2219892��
�

�
>__inference_c_layer_call_and_return_conditional_losses_2217144

inputs0
matmul_readvariableop_resource: -
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

: *
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:���������w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:��������� : : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:��������� 
 
_user_specified_nameinputs
�

�
>__inference_T_layer_call_and_return_conditional_losses_2217195

inputs0
matmul_readvariableop_resource: -
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

: *
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:���������w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:��������� : : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:��������� 
 
_user_specified_nameinputs
�^
�
B__inference_model_layer_call_and_return_conditional_losses_2217686

inputs
normalization_sub_y
normalization_sqrt_x

l2_2217587:@

l2_2217589:@

l3_2217592:@@

l3_2217594:@
	x_2217597:@@
	x_2217599:@
unknown:@@
	unknown_0:@
	unknown_1:@@
	unknown_2:@
	unknown_3:@@
	unknown_4:@
	unknown_5:@@
	unknown_6:@
	unknown_7:@@
	unknown_8:@
	unknown_9:@@

unknown_10:@

unknown_11:@@

unknown_12:@

unknown_13:@@

unknown_14:@

unknown_15:@ 

unknown_16: 

unknown_17:@ 

unknown_18: 

unknown_19:@ 

unknown_20: 

unknown_21:@ 

unknown_22: 
	c_2217662: 
	c_2217664:
phi_2217667: 
phi_2217669:
	p_2217672: 
	p_2217674:
	t_2217677: 
	t_2217679:
identity

identity_1

identity_2

identity_3��1/StatefulPartitionedCall�10/StatefulPartitionedCall�11/StatefulPartitionedCall�12/StatefulPartitionedCall�2/StatefulPartitionedCall�3/StatefulPartitionedCall�4/StatefulPartitionedCall�5/StatefulPartitionedCall�6/StatefulPartitionedCall�7/StatefulPartitionedCall�8/StatefulPartitionedCall�9/StatefulPartitionedCall�P/StatefulPartitionedCall�T/StatefulPartitionedCall�c/StatefulPartitionedCall�l2/StatefulPartitionedCall�l3/StatefulPartitionedCall�phi/StatefulPartitionedCall�x/StatefulPartitionedCallg
normalization/subSubinputsnormalization_sub_y*
T0*'
_output_shapes
:���������Y
normalization/SqrtSqrtnormalization_sqrt_x*
T0*
_output_shapes

:\
normalization/Maximum/yConst*
_output_shapes
: *
dtype0*
valueB
 *���3�
normalization/MaximumMaximumnormalization/Sqrt:y:0 normalization/Maximum/y:output:0*
T0*
_output_shapes

:�
normalization/truedivRealDivnormalization/sub:z:0normalization/Maximum:z:0*
T0*'
_output_shapes
:����������
l2/StatefulPartitionedCallStatefulPartitionedCallnormalization/truediv:z:0
l2_2217587
l2_2217589*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *H
fCRA
?__inference_l2_layer_call_and_return_conditional_losses_2216889�
l3/StatefulPartitionedCallStatefulPartitionedCall#l2/StatefulPartitionedCall:output:0
l3_2217592
l3_2217594*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *H
fCRA
?__inference_l3_layer_call_and_return_conditional_losses_2216906�
x/StatefulPartitionedCallStatefulPartitionedCall#l3/StatefulPartitionedCall:output:0	x_2217597	x_2217599*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_x_layer_call_and_return_conditional_losses_2216923�
10/StatefulPartitionedCallStatefulPartitionedCall"x/StatefulPartitionedCall:output:0unknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *H
fCRA
?__inference_10_layer_call_and_return_conditional_losses_2216940�
7/StatefulPartitionedCallStatefulPartitionedCall"x/StatefulPartitionedCall:output:0	unknown_1	unknown_2*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_7_layer_call_and_return_conditional_losses_2216957�
4/StatefulPartitionedCallStatefulPartitionedCall"x/StatefulPartitionedCall:output:0	unknown_3	unknown_4*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_4_layer_call_and_return_conditional_losses_2216974�
1/StatefulPartitionedCallStatefulPartitionedCall"x/StatefulPartitionedCall:output:0	unknown_5	unknown_6*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_1_layer_call_and_return_conditional_losses_2216991�
11/StatefulPartitionedCallStatefulPartitionedCall#10/StatefulPartitionedCall:output:0	unknown_7	unknown_8*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *H
fCRA
?__inference_11_layer_call_and_return_conditional_losses_2217008�
8/StatefulPartitionedCallStatefulPartitionedCall"7/StatefulPartitionedCall:output:0	unknown_9
unknown_10*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_8_layer_call_and_return_conditional_losses_2217025�
5/StatefulPartitionedCallStatefulPartitionedCall"4/StatefulPartitionedCall:output:0
unknown_11
unknown_12*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_5_layer_call_and_return_conditional_losses_2217042�
2/StatefulPartitionedCallStatefulPartitionedCall"1/StatefulPartitionedCall:output:0
unknown_13
unknown_14*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_2_layer_call_and_return_conditional_losses_2217059�
12/StatefulPartitionedCallStatefulPartitionedCall#11/StatefulPartitionedCall:output:0
unknown_15
unknown_16*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:��������� *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *H
fCRA
?__inference_12_layer_call_and_return_conditional_losses_2217076�
9/StatefulPartitionedCallStatefulPartitionedCall"8/StatefulPartitionedCall:output:0
unknown_17
unknown_18*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:��������� *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_9_layer_call_and_return_conditional_losses_2217093�
6/StatefulPartitionedCallStatefulPartitionedCall"5/StatefulPartitionedCall:output:0
unknown_19
unknown_20*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:��������� *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_6_layer_call_and_return_conditional_losses_2217110�
3/StatefulPartitionedCallStatefulPartitionedCall"2/StatefulPartitionedCall:output:0
unknown_21
unknown_22*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:��������� *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_3_layer_call_and_return_conditional_losses_2217127�
c/StatefulPartitionedCallStatefulPartitionedCall#12/StatefulPartitionedCall:output:0	c_2217662	c_2217664*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_c_layer_call_and_return_conditional_losses_2217144�
phi/StatefulPartitionedCallStatefulPartitionedCall"9/StatefulPartitionedCall:output:0phi_2217667phi_2217669*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *I
fDRB
@__inference_phi_layer_call_and_return_conditional_losses_2217161�
P/StatefulPartitionedCallStatefulPartitionedCall"6/StatefulPartitionedCall:output:0	p_2217672	p_2217674*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_P_layer_call_and_return_conditional_losses_2217178�
T/StatefulPartitionedCallStatefulPartitionedCall"3/StatefulPartitionedCall:output:0	t_2217677	t_2217679*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_T_layer_call_and_return_conditional_losses_2217195q
IdentityIdentity"T/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������s

Identity_1Identity"P/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������u

Identity_2Identity$phi/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������s

Identity_3Identity"c/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp^1/StatefulPartitionedCall^10/StatefulPartitionedCall^11/StatefulPartitionedCall^12/StatefulPartitionedCall^2/StatefulPartitionedCall^3/StatefulPartitionedCall^4/StatefulPartitionedCall^5/StatefulPartitionedCall^6/StatefulPartitionedCall^7/StatefulPartitionedCall^8/StatefulPartitionedCall^9/StatefulPartitionedCall^P/StatefulPartitionedCall^T/StatefulPartitionedCall^c/StatefulPartitionedCall^l2/StatefulPartitionedCall^l3/StatefulPartitionedCall^phi/StatefulPartitionedCall^x/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0"!

identity_3Identity_3:output:0*(
_construction_contextkEagerRuntime*�
_input_shapesu
s:���������::: : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : 26
1/StatefulPartitionedCall1/StatefulPartitionedCall28
10/StatefulPartitionedCall10/StatefulPartitionedCall28
11/StatefulPartitionedCall11/StatefulPartitionedCall28
12/StatefulPartitionedCall12/StatefulPartitionedCall26
2/StatefulPartitionedCall2/StatefulPartitionedCall26
3/StatefulPartitionedCall3/StatefulPartitionedCall26
4/StatefulPartitionedCall4/StatefulPartitionedCall26
5/StatefulPartitionedCall5/StatefulPartitionedCall26
6/StatefulPartitionedCall6/StatefulPartitionedCall26
7/StatefulPartitionedCall7/StatefulPartitionedCall26
8/StatefulPartitionedCall8/StatefulPartitionedCall26
9/StatefulPartitionedCall9/StatefulPartitionedCall26
P/StatefulPartitionedCallP/StatefulPartitionedCall26
T/StatefulPartitionedCallT/StatefulPartitionedCall26
c/StatefulPartitionedCallc/StatefulPartitionedCall28
l2/StatefulPartitionedCalll2/StatefulPartitionedCall28
l3/StatefulPartitionedCalll3/StatefulPartitionedCall2:
phi/StatefulPartitionedCallphi/StatefulPartitionedCall26
x/StatefulPartitionedCallx/StatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs:$ 

_output_shapes

::$ 

_output_shapes

:
�^
�
B__inference_model_layer_call_and_return_conditional_losses_2217205

inputs
normalization_sub_y
normalization_sqrt_x

l2_2216890:@

l2_2216892:@

l3_2216907:@@

l3_2216909:@
	x_2216924:@@
	x_2216926:@
unknown:@@
	unknown_0:@
	unknown_1:@@
	unknown_2:@
	unknown_3:@@
	unknown_4:@
	unknown_5:@@
	unknown_6:@
	unknown_7:@@
	unknown_8:@
	unknown_9:@@

unknown_10:@

unknown_11:@@

unknown_12:@

unknown_13:@@

unknown_14:@

unknown_15:@ 

unknown_16: 

unknown_17:@ 

unknown_18: 

unknown_19:@ 

unknown_20: 

unknown_21:@ 

unknown_22: 
	c_2217145: 
	c_2217147:
phi_2217162: 
phi_2217164:
	p_2217179: 
	p_2217181:
	t_2217196: 
	t_2217198:
identity

identity_1

identity_2

identity_3��1/StatefulPartitionedCall�10/StatefulPartitionedCall�11/StatefulPartitionedCall�12/StatefulPartitionedCall�2/StatefulPartitionedCall�3/StatefulPartitionedCall�4/StatefulPartitionedCall�5/StatefulPartitionedCall�6/StatefulPartitionedCall�7/StatefulPartitionedCall�8/StatefulPartitionedCall�9/StatefulPartitionedCall�P/StatefulPartitionedCall�T/StatefulPartitionedCall�c/StatefulPartitionedCall�l2/StatefulPartitionedCall�l3/StatefulPartitionedCall�phi/StatefulPartitionedCall�x/StatefulPartitionedCallg
normalization/subSubinputsnormalization_sub_y*
T0*'
_output_shapes
:���������Y
normalization/SqrtSqrtnormalization_sqrt_x*
T0*
_output_shapes

:\
normalization/Maximum/yConst*
_output_shapes
: *
dtype0*
valueB
 *���3�
normalization/MaximumMaximumnormalization/Sqrt:y:0 normalization/Maximum/y:output:0*
T0*
_output_shapes

:�
normalization/truedivRealDivnormalization/sub:z:0normalization/Maximum:z:0*
T0*'
_output_shapes
:����������
l2/StatefulPartitionedCallStatefulPartitionedCallnormalization/truediv:z:0
l2_2216890
l2_2216892*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *H
fCRA
?__inference_l2_layer_call_and_return_conditional_losses_2216889�
l3/StatefulPartitionedCallStatefulPartitionedCall#l2/StatefulPartitionedCall:output:0
l3_2216907
l3_2216909*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *H
fCRA
?__inference_l3_layer_call_and_return_conditional_losses_2216906�
x/StatefulPartitionedCallStatefulPartitionedCall#l3/StatefulPartitionedCall:output:0	x_2216924	x_2216926*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_x_layer_call_and_return_conditional_losses_2216923�
10/StatefulPartitionedCallStatefulPartitionedCall"x/StatefulPartitionedCall:output:0unknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *H
fCRA
?__inference_10_layer_call_and_return_conditional_losses_2216940�
7/StatefulPartitionedCallStatefulPartitionedCall"x/StatefulPartitionedCall:output:0	unknown_1	unknown_2*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_7_layer_call_and_return_conditional_losses_2216957�
4/StatefulPartitionedCallStatefulPartitionedCall"x/StatefulPartitionedCall:output:0	unknown_3	unknown_4*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_4_layer_call_and_return_conditional_losses_2216974�
1/StatefulPartitionedCallStatefulPartitionedCall"x/StatefulPartitionedCall:output:0	unknown_5	unknown_6*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_1_layer_call_and_return_conditional_losses_2216991�
11/StatefulPartitionedCallStatefulPartitionedCall#10/StatefulPartitionedCall:output:0	unknown_7	unknown_8*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *H
fCRA
?__inference_11_layer_call_and_return_conditional_losses_2217008�
8/StatefulPartitionedCallStatefulPartitionedCall"7/StatefulPartitionedCall:output:0	unknown_9
unknown_10*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_8_layer_call_and_return_conditional_losses_2217025�
5/StatefulPartitionedCallStatefulPartitionedCall"4/StatefulPartitionedCall:output:0
unknown_11
unknown_12*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_5_layer_call_and_return_conditional_losses_2217042�
2/StatefulPartitionedCallStatefulPartitionedCall"1/StatefulPartitionedCall:output:0
unknown_13
unknown_14*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_2_layer_call_and_return_conditional_losses_2217059�
12/StatefulPartitionedCallStatefulPartitionedCall#11/StatefulPartitionedCall:output:0
unknown_15
unknown_16*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:��������� *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *H
fCRA
?__inference_12_layer_call_and_return_conditional_losses_2217076�
9/StatefulPartitionedCallStatefulPartitionedCall"8/StatefulPartitionedCall:output:0
unknown_17
unknown_18*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:��������� *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_9_layer_call_and_return_conditional_losses_2217093�
6/StatefulPartitionedCallStatefulPartitionedCall"5/StatefulPartitionedCall:output:0
unknown_19
unknown_20*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:��������� *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_6_layer_call_and_return_conditional_losses_2217110�
3/StatefulPartitionedCallStatefulPartitionedCall"2/StatefulPartitionedCall:output:0
unknown_21
unknown_22*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:��������� *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_3_layer_call_and_return_conditional_losses_2217127�
c/StatefulPartitionedCallStatefulPartitionedCall#12/StatefulPartitionedCall:output:0	c_2217145	c_2217147*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_c_layer_call_and_return_conditional_losses_2217144�
phi/StatefulPartitionedCallStatefulPartitionedCall"9/StatefulPartitionedCall:output:0phi_2217162phi_2217164*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *I
fDRB
@__inference_phi_layer_call_and_return_conditional_losses_2217161�
P/StatefulPartitionedCallStatefulPartitionedCall"6/StatefulPartitionedCall:output:0	p_2217179	p_2217181*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_P_layer_call_and_return_conditional_losses_2217178�
T/StatefulPartitionedCallStatefulPartitionedCall"3/StatefulPartitionedCall:output:0	t_2217196	t_2217198*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_T_layer_call_and_return_conditional_losses_2217195q
IdentityIdentity"T/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������s

Identity_1Identity"P/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������u

Identity_2Identity$phi/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������s

Identity_3Identity"c/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp^1/StatefulPartitionedCall^10/StatefulPartitionedCall^11/StatefulPartitionedCall^12/StatefulPartitionedCall^2/StatefulPartitionedCall^3/StatefulPartitionedCall^4/StatefulPartitionedCall^5/StatefulPartitionedCall^6/StatefulPartitionedCall^7/StatefulPartitionedCall^8/StatefulPartitionedCall^9/StatefulPartitionedCall^P/StatefulPartitionedCall^T/StatefulPartitionedCall^c/StatefulPartitionedCall^l2/StatefulPartitionedCall^l3/StatefulPartitionedCall^phi/StatefulPartitionedCall^x/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0"!

identity_3Identity_3:output:0*(
_construction_contextkEagerRuntime*�
_input_shapesu
s:���������::: : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : 26
1/StatefulPartitionedCall1/StatefulPartitionedCall28
10/StatefulPartitionedCall10/StatefulPartitionedCall28
11/StatefulPartitionedCall11/StatefulPartitionedCall28
12/StatefulPartitionedCall12/StatefulPartitionedCall26
2/StatefulPartitionedCall2/StatefulPartitionedCall26
3/StatefulPartitionedCall3/StatefulPartitionedCall26
4/StatefulPartitionedCall4/StatefulPartitionedCall26
5/StatefulPartitionedCall5/StatefulPartitionedCall26
6/StatefulPartitionedCall6/StatefulPartitionedCall26
7/StatefulPartitionedCall7/StatefulPartitionedCall26
8/StatefulPartitionedCall8/StatefulPartitionedCall26
9/StatefulPartitionedCall9/StatefulPartitionedCall26
P/StatefulPartitionedCallP/StatefulPartitionedCall26
T/StatefulPartitionedCallT/StatefulPartitionedCall26
c/StatefulPartitionedCallc/StatefulPartitionedCall28
l2/StatefulPartitionedCalll2/StatefulPartitionedCall28
l3/StatefulPartitionedCalll3/StatefulPartitionedCall2:
phi/StatefulPartitionedCallphi/StatefulPartitionedCall26
x/StatefulPartitionedCallx/StatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs:$ 

_output_shapes

::$ 

_output_shapes

:
�
�
$__inference_10_layer_call_fn_2218829

inputs
unknown:@@
	unknown_0:@
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *H
fCRA
?__inference_10_layer_call_and_return_conditional_losses_2216940o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������@`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�

�
?__inference_10_layer_call_and_return_conditional_losses_2218840

inputs0
matmul_readvariableop_resource:@@-
biasadd_readvariableop_resource:@
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@@*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������@a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:���������@w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�
�	
%__inference_signature_wrapper_2218179	
input
unknown
	unknown_0
	unknown_1:@
	unknown_2:@
	unknown_3:@@
	unknown_4:@
	unknown_5:@@
	unknown_6:@
	unknown_7:@@
	unknown_8:@
	unknown_9:@@

unknown_10:@

unknown_11:@@

unknown_12:@

unknown_13:@@

unknown_14:@

unknown_15:@@

unknown_16:@

unknown_17:@@

unknown_18:@

unknown_19:@@

unknown_20:@

unknown_21:@@

unknown_22:@

unknown_23:@ 

unknown_24: 

unknown_25:@ 

unknown_26: 

unknown_27:@ 

unknown_28: 

unknown_29:@ 

unknown_30: 

unknown_31: 

unknown_32:

unknown_33: 

unknown_34:

unknown_35: 

unknown_36:

unknown_37: 

unknown_38:
identity

identity_1

identity_2

identity_3��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17
unknown_18
unknown_19
unknown_20
unknown_21
unknown_22
unknown_23
unknown_24
unknown_25
unknown_26
unknown_27
unknown_28
unknown_29
unknown_30
unknown_31
unknown_32
unknown_33
unknown_34
unknown_35
unknown_36
unknown_37
unknown_38*4
Tin-
+2)*
Tout
2*
_collective_manager_ids
 *`
_output_shapesN
L:���������:���������:���������:���������*H
_read_only_resource_inputs*
(&	
 !"#$%&'(*-
config_proto

CPU

GPU 2J 8� *+
f&R$
"__inference__wrapped_model_2216864o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������q

Identity_1Identity StatefulPartitionedCall:output:1^NoOp*
T0*'
_output_shapes
:���������q

Identity_2Identity StatefulPartitionedCall:output:2^NoOp*
T0*'
_output_shapes
:���������q

Identity_3Identity StatefulPartitionedCall:output:3^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0"!

identity_3Identity_3:output:0*(
_construction_contextkEagerRuntime*�
_input_shapesu
s:���������::: : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:N J
'
_output_shapes
:���������

_user_specified_nameinput:$ 

_output_shapes

::$ 

_output_shapes

:
�

�
>__inference_5_layer_call_and_return_conditional_losses_2217042

inputs0
matmul_readvariableop_resource:@@-
biasadd_readvariableop_resource:@
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@@*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������@a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:���������@w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�

�
@__inference_phi_layer_call_and_return_conditional_losses_2217161

inputs0
matmul_readvariableop_resource: -
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

: *
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:���������w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:��������� : : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:��������� 
 
_user_specified_nameinputs
�

�
>__inference_T_layer_call_and_return_conditional_losses_2219020

inputs0
matmul_readvariableop_resource: -
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

: *
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:���������w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:��������� : : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:��������� 
 
_user_specified_nameinputs
�

�
>__inference_4_layer_call_and_return_conditional_losses_2216974

inputs0
matmul_readvariableop_resource:@@-
biasadd_readvariableop_resource:@
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@@*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������@a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:���������@w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�

�
>__inference_7_layer_call_and_return_conditional_losses_2216957

inputs0
matmul_readvariableop_resource:@@-
biasadd_readvariableop_resource:@
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@@*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������@a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:���������@w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�
�
$__inference_11_layer_call_fn_2218909

inputs
unknown:@@
	unknown_0:@
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *H
fCRA
?__inference_11_layer_call_and_return_conditional_losses_2217008o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������@`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�

�
>__inference_4_layer_call_and_return_conditional_losses_2218800

inputs0
matmul_readvariableop_resource:@@-
biasadd_readvariableop_resource:@
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@@*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������@a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:���������@w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
ߧ
�
"__inference__wrapped_model_2216864	
input
model_normalization_sub_y
model_normalization_sqrt_x9
'model_l2_matmul_readvariableop_resource:@6
(model_l2_biasadd_readvariableop_resource:@9
'model_l3_matmul_readvariableop_resource:@@6
(model_l3_biasadd_readvariableop_resource:@8
&model_x_matmul_readvariableop_resource:@@5
'model_x_biasadd_readvariableop_resource:@9
'model_10_matmul_readvariableop_resource:@@6
(model_10_biasadd_readvariableop_resource:@8
&model_7_matmul_readvariableop_resource:@@5
'model_7_biasadd_readvariableop_resource:@8
&model_4_matmul_readvariableop_resource:@@5
'model_4_biasadd_readvariableop_resource:@8
&model_1_matmul_readvariableop_resource:@@5
'model_1_biasadd_readvariableop_resource:@9
'model_11_matmul_readvariableop_resource:@@6
(model_11_biasadd_readvariableop_resource:@8
&model_8_matmul_readvariableop_resource:@@5
'model_8_biasadd_readvariableop_resource:@8
&model_5_matmul_readvariableop_resource:@@5
'model_5_biasadd_readvariableop_resource:@8
&model_2_matmul_readvariableop_resource:@@5
'model_2_biasadd_readvariableop_resource:@9
'model_12_matmul_readvariableop_resource:@ 6
(model_12_biasadd_readvariableop_resource: 8
&model_9_matmul_readvariableop_resource:@ 5
'model_9_biasadd_readvariableop_resource: 8
&model_6_matmul_readvariableop_resource:@ 5
'model_6_biasadd_readvariableop_resource: 8
&model_3_matmul_readvariableop_resource:@ 5
'model_3_biasadd_readvariableop_resource: 8
&model_c_matmul_readvariableop_resource: 5
'model_c_biasadd_readvariableop_resource::
(model_phi_matmul_readvariableop_resource: 7
)model_phi_biasadd_readvariableop_resource:8
&model_p_matmul_readvariableop_resource: 5
'model_p_biasadd_readvariableop_resource:8
&model_t_matmul_readvariableop_resource: 5
'model_t_biasadd_readvariableop_resource:
identity

identity_1

identity_2

identity_3��model/1/BiasAdd/ReadVariableOp�model/1/MatMul/ReadVariableOp�model/10/BiasAdd/ReadVariableOp�model/10/MatMul/ReadVariableOp�model/11/BiasAdd/ReadVariableOp�model/11/MatMul/ReadVariableOp�model/12/BiasAdd/ReadVariableOp�model/12/MatMul/ReadVariableOp�model/2/BiasAdd/ReadVariableOp�model/2/MatMul/ReadVariableOp�model/3/BiasAdd/ReadVariableOp�model/3/MatMul/ReadVariableOp�model/4/BiasAdd/ReadVariableOp�model/4/MatMul/ReadVariableOp�model/5/BiasAdd/ReadVariableOp�model/5/MatMul/ReadVariableOp�model/6/BiasAdd/ReadVariableOp�model/6/MatMul/ReadVariableOp�model/7/BiasAdd/ReadVariableOp�model/7/MatMul/ReadVariableOp�model/8/BiasAdd/ReadVariableOp�model/8/MatMul/ReadVariableOp�model/9/BiasAdd/ReadVariableOp�model/9/MatMul/ReadVariableOp�model/P/BiasAdd/ReadVariableOp�model/P/MatMul/ReadVariableOp�model/T/BiasAdd/ReadVariableOp�model/T/MatMul/ReadVariableOp�model/c/BiasAdd/ReadVariableOp�model/c/MatMul/ReadVariableOp�model/l2/BiasAdd/ReadVariableOp�model/l2/MatMul/ReadVariableOp�model/l3/BiasAdd/ReadVariableOp�model/l3/MatMul/ReadVariableOp� model/phi/BiasAdd/ReadVariableOp�model/phi/MatMul/ReadVariableOp�model/x/BiasAdd/ReadVariableOp�model/x/MatMul/ReadVariableOpr
model/normalization/subSubinputmodel_normalization_sub_y*
T0*'
_output_shapes
:���������e
model/normalization/SqrtSqrtmodel_normalization_sqrt_x*
T0*
_output_shapes

:b
model/normalization/Maximum/yConst*
_output_shapes
: *
dtype0*
valueB
 *���3�
model/normalization/MaximumMaximummodel/normalization/Sqrt:y:0&model/normalization/Maximum/y:output:0*
T0*
_output_shapes

:�
model/normalization/truedivRealDivmodel/normalization/sub:z:0model/normalization/Maximum:z:0*
T0*'
_output_shapes
:����������
model/l2/MatMul/ReadVariableOpReadVariableOp'model_l2_matmul_readvariableop_resource*
_output_shapes

:@*
dtype0�
model/l2/MatMulMatMulmodel/normalization/truediv:z:0&model/l2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@�
model/l2/BiasAdd/ReadVariableOpReadVariableOp(model_l2_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype0�
model/l2/BiasAddBiasAddmodel/l2/MatMul:product:0'model/l2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@b
model/l2/ReluRelumodel/l2/BiasAdd:output:0*
T0*'
_output_shapes
:���������@�
model/l3/MatMul/ReadVariableOpReadVariableOp'model_l3_matmul_readvariableop_resource*
_output_shapes

:@@*
dtype0�
model/l3/MatMulMatMulmodel/l2/Relu:activations:0&model/l3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@�
model/l3/BiasAdd/ReadVariableOpReadVariableOp(model_l3_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype0�
model/l3/BiasAddBiasAddmodel/l3/MatMul:product:0'model/l3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@b
model/l3/ReluRelumodel/l3/BiasAdd:output:0*
T0*'
_output_shapes
:���������@�
model/x/MatMul/ReadVariableOpReadVariableOp&model_x_matmul_readvariableop_resource*
_output_shapes

:@@*
dtype0�
model/x/MatMulMatMulmodel/l3/Relu:activations:0%model/x/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@�
model/x/BiasAdd/ReadVariableOpReadVariableOp'model_x_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype0�
model/x/BiasAddBiasAddmodel/x/MatMul:product:0&model/x/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@`
model/x/ReluRelumodel/x/BiasAdd:output:0*
T0*'
_output_shapes
:���������@�
model/10/MatMul/ReadVariableOpReadVariableOp'model_10_matmul_readvariableop_resource*
_output_shapes

:@@*
dtype0�
model/10/MatMulMatMulmodel/x/Relu:activations:0&model/10/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@�
model/10/BiasAdd/ReadVariableOpReadVariableOp(model_10_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype0�
model/10/BiasAddBiasAddmodel/10/MatMul:product:0'model/10/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@b
model/10/ReluRelumodel/10/BiasAdd:output:0*
T0*'
_output_shapes
:���������@�
model/7/MatMul/ReadVariableOpReadVariableOp&model_7_matmul_readvariableop_resource*
_output_shapes

:@@*
dtype0�
model/7/MatMulMatMulmodel/x/Relu:activations:0%model/7/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@�
model/7/BiasAdd/ReadVariableOpReadVariableOp'model_7_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype0�
model/7/BiasAddBiasAddmodel/7/MatMul:product:0&model/7/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@`
model/7/ReluRelumodel/7/BiasAdd:output:0*
T0*'
_output_shapes
:���������@�
model/4/MatMul/ReadVariableOpReadVariableOp&model_4_matmul_readvariableop_resource*
_output_shapes

:@@*
dtype0�
model/4/MatMulMatMulmodel/x/Relu:activations:0%model/4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@�
model/4/BiasAdd/ReadVariableOpReadVariableOp'model_4_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype0�
model/4/BiasAddBiasAddmodel/4/MatMul:product:0&model/4/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@`
model/4/ReluRelumodel/4/BiasAdd:output:0*
T0*'
_output_shapes
:���������@�
model/1/MatMul/ReadVariableOpReadVariableOp&model_1_matmul_readvariableop_resource*
_output_shapes

:@@*
dtype0�
model/1/MatMulMatMulmodel/x/Relu:activations:0%model/1/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@�
model/1/BiasAdd/ReadVariableOpReadVariableOp'model_1_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype0�
model/1/BiasAddBiasAddmodel/1/MatMul:product:0&model/1/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@`
model/1/ReluRelumodel/1/BiasAdd:output:0*
T0*'
_output_shapes
:���������@�
model/11/MatMul/ReadVariableOpReadVariableOp'model_11_matmul_readvariableop_resource*
_output_shapes

:@@*
dtype0�
model/11/MatMulMatMulmodel/10/Relu:activations:0&model/11/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@�
model/11/BiasAdd/ReadVariableOpReadVariableOp(model_11_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype0�
model/11/BiasAddBiasAddmodel/11/MatMul:product:0'model/11/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@b
model/11/ReluRelumodel/11/BiasAdd:output:0*
T0*'
_output_shapes
:���������@�
model/8/MatMul/ReadVariableOpReadVariableOp&model_8_matmul_readvariableop_resource*
_output_shapes

:@@*
dtype0�
model/8/MatMulMatMulmodel/7/Relu:activations:0%model/8/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@�
model/8/BiasAdd/ReadVariableOpReadVariableOp'model_8_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype0�
model/8/BiasAddBiasAddmodel/8/MatMul:product:0&model/8/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@`
model/8/ReluRelumodel/8/BiasAdd:output:0*
T0*'
_output_shapes
:���������@�
model/5/MatMul/ReadVariableOpReadVariableOp&model_5_matmul_readvariableop_resource*
_output_shapes

:@@*
dtype0�
model/5/MatMulMatMulmodel/4/Relu:activations:0%model/5/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@�
model/5/BiasAdd/ReadVariableOpReadVariableOp'model_5_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype0�
model/5/BiasAddBiasAddmodel/5/MatMul:product:0&model/5/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@`
model/5/ReluRelumodel/5/BiasAdd:output:0*
T0*'
_output_shapes
:���������@�
model/2/MatMul/ReadVariableOpReadVariableOp&model_2_matmul_readvariableop_resource*
_output_shapes

:@@*
dtype0�
model/2/MatMulMatMulmodel/1/Relu:activations:0%model/2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@�
model/2/BiasAdd/ReadVariableOpReadVariableOp'model_2_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype0�
model/2/BiasAddBiasAddmodel/2/MatMul:product:0&model/2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@`
model/2/ReluRelumodel/2/BiasAdd:output:0*
T0*'
_output_shapes
:���������@�
model/12/MatMul/ReadVariableOpReadVariableOp'model_12_matmul_readvariableop_resource*
_output_shapes

:@ *
dtype0�
model/12/MatMulMatMulmodel/11/Relu:activations:0&model/12/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� �
model/12/BiasAdd/ReadVariableOpReadVariableOp(model_12_biasadd_readvariableop_resource*
_output_shapes
: *
dtype0�
model/12/BiasAddBiasAddmodel/12/MatMul:product:0'model/12/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� b
model/12/ReluRelumodel/12/BiasAdd:output:0*
T0*'
_output_shapes
:��������� �
model/9/MatMul/ReadVariableOpReadVariableOp&model_9_matmul_readvariableop_resource*
_output_shapes

:@ *
dtype0�
model/9/MatMulMatMulmodel/8/Relu:activations:0%model/9/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� �
model/9/BiasAdd/ReadVariableOpReadVariableOp'model_9_biasadd_readvariableop_resource*
_output_shapes
: *
dtype0�
model/9/BiasAddBiasAddmodel/9/MatMul:product:0&model/9/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� `
model/9/ReluRelumodel/9/BiasAdd:output:0*
T0*'
_output_shapes
:��������� �
model/6/MatMul/ReadVariableOpReadVariableOp&model_6_matmul_readvariableop_resource*
_output_shapes

:@ *
dtype0�
model/6/MatMulMatMulmodel/5/Relu:activations:0%model/6/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� �
model/6/BiasAdd/ReadVariableOpReadVariableOp'model_6_biasadd_readvariableop_resource*
_output_shapes
: *
dtype0�
model/6/BiasAddBiasAddmodel/6/MatMul:product:0&model/6/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� `
model/6/ReluRelumodel/6/BiasAdd:output:0*
T0*'
_output_shapes
:��������� �
model/3/MatMul/ReadVariableOpReadVariableOp&model_3_matmul_readvariableop_resource*
_output_shapes

:@ *
dtype0�
model/3/MatMulMatMulmodel/2/Relu:activations:0%model/3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� �
model/3/BiasAdd/ReadVariableOpReadVariableOp'model_3_biasadd_readvariableop_resource*
_output_shapes
: *
dtype0�
model/3/BiasAddBiasAddmodel/3/MatMul:product:0&model/3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� `
model/3/ReluRelumodel/3/BiasAdd:output:0*
T0*'
_output_shapes
:��������� �
model/c/MatMul/ReadVariableOpReadVariableOp&model_c_matmul_readvariableop_resource*
_output_shapes

: *
dtype0�
model/c/MatMulMatMulmodel/12/Relu:activations:0%model/c/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
model/c/BiasAdd/ReadVariableOpReadVariableOp'model_c_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
model/c/BiasAddBiasAddmodel/c/MatMul:product:0&model/c/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������`
model/c/ReluRelumodel/c/BiasAdd:output:0*
T0*'
_output_shapes
:����������
model/phi/MatMul/ReadVariableOpReadVariableOp(model_phi_matmul_readvariableop_resource*
_output_shapes

: *
dtype0�
model/phi/MatMulMatMulmodel/9/Relu:activations:0'model/phi/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
 model/phi/BiasAdd/ReadVariableOpReadVariableOp)model_phi_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
model/phi/BiasAddBiasAddmodel/phi/MatMul:product:0(model/phi/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������d
model/phi/ReluRelumodel/phi/BiasAdd:output:0*
T0*'
_output_shapes
:����������
model/P/MatMul/ReadVariableOpReadVariableOp&model_p_matmul_readvariableop_resource*
_output_shapes

: *
dtype0�
model/P/MatMulMatMulmodel/6/Relu:activations:0%model/P/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
model/P/BiasAdd/ReadVariableOpReadVariableOp'model_p_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
model/P/BiasAddBiasAddmodel/P/MatMul:product:0&model/P/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������`
model/P/ReluRelumodel/P/BiasAdd:output:0*
T0*'
_output_shapes
:����������
model/T/MatMul/ReadVariableOpReadVariableOp&model_t_matmul_readvariableop_resource*
_output_shapes

: *
dtype0�
model/T/MatMulMatMulmodel/3/Relu:activations:0%model/T/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:����������
model/T/BiasAdd/ReadVariableOpReadVariableOp'model_t_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
model/T/BiasAddBiasAddmodel/T/MatMul:product:0&model/T/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������`
model/T/ReluRelumodel/T/BiasAdd:output:0*
T0*'
_output_shapes
:���������i
IdentityIdentitymodel/P/Relu:activations:0^NoOp*
T0*'
_output_shapes
:���������k

Identity_1Identitymodel/T/Relu:activations:0^NoOp*
T0*'
_output_shapes
:���������k

Identity_2Identitymodel/c/Relu:activations:0^NoOp*
T0*'
_output_shapes
:���������m

Identity_3Identitymodel/phi/Relu:activations:0^NoOp*
T0*'
_output_shapes
:����������

NoOpNoOp^model/1/BiasAdd/ReadVariableOp^model/1/MatMul/ReadVariableOp ^model/10/BiasAdd/ReadVariableOp^model/10/MatMul/ReadVariableOp ^model/11/BiasAdd/ReadVariableOp^model/11/MatMul/ReadVariableOp ^model/12/BiasAdd/ReadVariableOp^model/12/MatMul/ReadVariableOp^model/2/BiasAdd/ReadVariableOp^model/2/MatMul/ReadVariableOp^model/3/BiasAdd/ReadVariableOp^model/3/MatMul/ReadVariableOp^model/4/BiasAdd/ReadVariableOp^model/4/MatMul/ReadVariableOp^model/5/BiasAdd/ReadVariableOp^model/5/MatMul/ReadVariableOp^model/6/BiasAdd/ReadVariableOp^model/6/MatMul/ReadVariableOp^model/7/BiasAdd/ReadVariableOp^model/7/MatMul/ReadVariableOp^model/8/BiasAdd/ReadVariableOp^model/8/MatMul/ReadVariableOp^model/9/BiasAdd/ReadVariableOp^model/9/MatMul/ReadVariableOp^model/P/BiasAdd/ReadVariableOp^model/P/MatMul/ReadVariableOp^model/T/BiasAdd/ReadVariableOp^model/T/MatMul/ReadVariableOp^model/c/BiasAdd/ReadVariableOp^model/c/MatMul/ReadVariableOp ^model/l2/BiasAdd/ReadVariableOp^model/l2/MatMul/ReadVariableOp ^model/l3/BiasAdd/ReadVariableOp^model/l3/MatMul/ReadVariableOp!^model/phi/BiasAdd/ReadVariableOp ^model/phi/MatMul/ReadVariableOp^model/x/BiasAdd/ReadVariableOp^model/x/MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0"!

identity_3Identity_3:output:0*(
_construction_contextkEagerRuntime*�
_input_shapesu
s:���������::: : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : 2@
model/1/BiasAdd/ReadVariableOpmodel/1/BiasAdd/ReadVariableOp2>
model/1/MatMul/ReadVariableOpmodel/1/MatMul/ReadVariableOp2B
model/10/BiasAdd/ReadVariableOpmodel/10/BiasAdd/ReadVariableOp2@
model/10/MatMul/ReadVariableOpmodel/10/MatMul/ReadVariableOp2B
model/11/BiasAdd/ReadVariableOpmodel/11/BiasAdd/ReadVariableOp2@
model/11/MatMul/ReadVariableOpmodel/11/MatMul/ReadVariableOp2B
model/12/BiasAdd/ReadVariableOpmodel/12/BiasAdd/ReadVariableOp2@
model/12/MatMul/ReadVariableOpmodel/12/MatMul/ReadVariableOp2@
model/2/BiasAdd/ReadVariableOpmodel/2/BiasAdd/ReadVariableOp2>
model/2/MatMul/ReadVariableOpmodel/2/MatMul/ReadVariableOp2@
model/3/BiasAdd/ReadVariableOpmodel/3/BiasAdd/ReadVariableOp2>
model/3/MatMul/ReadVariableOpmodel/3/MatMul/ReadVariableOp2@
model/4/BiasAdd/ReadVariableOpmodel/4/BiasAdd/ReadVariableOp2>
model/4/MatMul/ReadVariableOpmodel/4/MatMul/ReadVariableOp2@
model/5/BiasAdd/ReadVariableOpmodel/5/BiasAdd/ReadVariableOp2>
model/5/MatMul/ReadVariableOpmodel/5/MatMul/ReadVariableOp2@
model/6/BiasAdd/ReadVariableOpmodel/6/BiasAdd/ReadVariableOp2>
model/6/MatMul/ReadVariableOpmodel/6/MatMul/ReadVariableOp2@
model/7/BiasAdd/ReadVariableOpmodel/7/BiasAdd/ReadVariableOp2>
model/7/MatMul/ReadVariableOpmodel/7/MatMul/ReadVariableOp2@
model/8/BiasAdd/ReadVariableOpmodel/8/BiasAdd/ReadVariableOp2>
model/8/MatMul/ReadVariableOpmodel/8/MatMul/ReadVariableOp2@
model/9/BiasAdd/ReadVariableOpmodel/9/BiasAdd/ReadVariableOp2>
model/9/MatMul/ReadVariableOpmodel/9/MatMul/ReadVariableOp2@
model/P/BiasAdd/ReadVariableOpmodel/P/BiasAdd/ReadVariableOp2>
model/P/MatMul/ReadVariableOpmodel/P/MatMul/ReadVariableOp2@
model/T/BiasAdd/ReadVariableOpmodel/T/BiasAdd/ReadVariableOp2>
model/T/MatMul/ReadVariableOpmodel/T/MatMul/ReadVariableOp2@
model/c/BiasAdd/ReadVariableOpmodel/c/BiasAdd/ReadVariableOp2>
model/c/MatMul/ReadVariableOpmodel/c/MatMul/ReadVariableOp2B
model/l2/BiasAdd/ReadVariableOpmodel/l2/BiasAdd/ReadVariableOp2@
model/l2/MatMul/ReadVariableOpmodel/l2/MatMul/ReadVariableOp2B
model/l3/BiasAdd/ReadVariableOpmodel/l3/BiasAdd/ReadVariableOp2@
model/l3/MatMul/ReadVariableOpmodel/l3/MatMul/ReadVariableOp2D
 model/phi/BiasAdd/ReadVariableOp model/phi/BiasAdd/ReadVariableOp2B
model/phi/MatMul/ReadVariableOpmodel/phi/MatMul/ReadVariableOp2@
model/x/BiasAdd/ReadVariableOpmodel/x/BiasAdd/ReadVariableOp2>
model/x/MatMul/ReadVariableOpmodel/x/MatMul/ReadVariableOp:N J
'
_output_shapes
:���������

_user_specified_nameinput:$ 

_output_shapes

::$ 

_output_shapes

:
�

�
?__inference_l3_layer_call_and_return_conditional_losses_2216906

inputs0
matmul_readvariableop_resource:@@-
biasadd_readvariableop_resource:@
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@@*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������@a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:���������@w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
��
�
B__inference_model_layer_call_and_return_conditional_losses_2218553

inputs
normalization_sub_y
normalization_sqrt_x3
!l2_matmul_readvariableop_resource:@0
"l2_biasadd_readvariableop_resource:@3
!l3_matmul_readvariableop_resource:@@0
"l3_biasadd_readvariableop_resource:@2
 x_matmul_readvariableop_resource:@@/
!x_biasadd_readvariableop_resource:@0
matmul_readvariableop_resource:@@-
biasadd_readvariableop_resource:@2
 matmul_readvariableop_resource_0:@@/
!biasadd_readvariableop_resource_0:@2
 matmul_readvariableop_resource_1:@@/
!biasadd_readvariableop_resource_1:@2
 matmul_readvariableop_resource_2:@@/
!biasadd_readvariableop_resource_2:@2
 matmul_readvariableop_resource_3:@@/
!biasadd_readvariableop_resource_3:@2
 matmul_readvariableop_resource_4:@@/
!biasadd_readvariableop_resource_4:@2
 matmul_readvariableop_resource_5:@@/
!biasadd_readvariableop_resource_5:@2
 matmul_readvariableop_resource_6:@@/
!biasadd_readvariableop_resource_6:@2
 matmul_readvariableop_resource_7:@ /
!biasadd_readvariableop_resource_7: 2
 matmul_readvariableop_resource_8:@ /
!biasadd_readvariableop_resource_8: 2
 matmul_readvariableop_resource_9:@ /
!biasadd_readvariableop_resource_9: 3
!matmul_readvariableop_resource_10:@ 0
"biasadd_readvariableop_resource_10: 2
 c_matmul_readvariableop_resource: /
!c_biasadd_readvariableop_resource:4
"phi_matmul_readvariableop_resource: 1
#phi_biasadd_readvariableop_resource:2
 p_matmul_readvariableop_resource: /
!p_biasadd_readvariableop_resource:2
 t_matmul_readvariableop_resource: /
!t_biasadd_readvariableop_resource:
identity

identity_1

identity_2

identity_3��1/BiasAdd/ReadVariableOp�1/MatMul/ReadVariableOp�10/BiasAdd/ReadVariableOp�10/MatMul/ReadVariableOp�11/BiasAdd/ReadVariableOp�11/MatMul/ReadVariableOp�12/BiasAdd/ReadVariableOp�12/MatMul/ReadVariableOp�2/BiasAdd/ReadVariableOp�2/MatMul/ReadVariableOp�3/BiasAdd/ReadVariableOp�3/MatMul/ReadVariableOp�4/BiasAdd/ReadVariableOp�4/MatMul/ReadVariableOp�5/BiasAdd/ReadVariableOp�5/MatMul/ReadVariableOp�6/BiasAdd/ReadVariableOp�6/MatMul/ReadVariableOp�7/BiasAdd/ReadVariableOp�7/MatMul/ReadVariableOp�8/BiasAdd/ReadVariableOp�8/MatMul/ReadVariableOp�9/BiasAdd/ReadVariableOp�9/MatMul/ReadVariableOp�P/BiasAdd/ReadVariableOp�P/MatMul/ReadVariableOp�T/BiasAdd/ReadVariableOp�T/MatMul/ReadVariableOp�c/BiasAdd/ReadVariableOp�c/MatMul/ReadVariableOp�l2/BiasAdd/ReadVariableOp�l2/MatMul/ReadVariableOp�l3/BiasAdd/ReadVariableOp�l3/MatMul/ReadVariableOp�phi/BiasAdd/ReadVariableOp�phi/MatMul/ReadVariableOp�x/BiasAdd/ReadVariableOp�x/MatMul/ReadVariableOpg
normalization/subSubinputsnormalization_sub_y*
T0*'
_output_shapes
:���������Y
normalization/SqrtSqrtnormalization_sqrt_x*
T0*
_output_shapes

:\
normalization/Maximum/yConst*
_output_shapes
: *
dtype0*
valueB
 *���3�
normalization/MaximumMaximumnormalization/Sqrt:y:0 normalization/Maximum/y:output:0*
T0*
_output_shapes

:�
normalization/truedivRealDivnormalization/sub:z:0normalization/Maximum:z:0*
T0*'
_output_shapes
:���������z
l2/MatMul/ReadVariableOpReadVariableOp!l2_matmul_readvariableop_resource*
_output_shapes

:@*
dtype0�
	l2/MatMulMatMulnormalization/truediv:z:0 l2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@x
l2/BiasAdd/ReadVariableOpReadVariableOp"l2_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype0

l2/BiasAddBiasAddl2/MatMul:product:0!l2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@V
l2/ReluRelul2/BiasAdd:output:0*
T0*'
_output_shapes
:���������@z
l3/MatMul/ReadVariableOpReadVariableOp!l3_matmul_readvariableop_resource*
_output_shapes

:@@*
dtype0~
	l3/MatMulMatMull2/Relu:activations:0 l3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@x
l3/BiasAdd/ReadVariableOpReadVariableOp"l3_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype0

l3/BiasAddBiasAddl3/MatMul:product:0!l3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@V
l3/ReluRelul3/BiasAdd:output:0*
T0*'
_output_shapes
:���������@x
x/MatMul/ReadVariableOpReadVariableOp x_matmul_readvariableop_resource*
_output_shapes

:@@*
dtype0|
x/MatMulMatMull3/Relu:activations:0x/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@v
x/BiasAdd/ReadVariableOpReadVariableOp!x_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype0|
	x/BiasAddBiasAddx/MatMul:product:0 x/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@T
x/ReluRelux/BiasAdd:output:0*
T0*'
_output_shapes
:���������@w
10/MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@@*
dtype0}
	10/MatMulMatMulx/Relu:activations:0 10/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@u
10/BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype0

10/BiasAddBiasAdd10/MatMul:product:0!10/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@V
10/ReluRelu10/BiasAdd:output:0*
T0*'
_output_shapes
:���������@x
7/MatMul/ReadVariableOpReadVariableOp matmul_readvariableop_resource_0*
_output_shapes

:@@*
dtype0{
7/MatMulMatMulx/Relu:activations:07/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@v
7/BiasAdd/ReadVariableOpReadVariableOp!biasadd_readvariableop_resource_0*
_output_shapes
:@*
dtype0|
	7/BiasAddBiasAdd7/MatMul:product:0 7/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@T
7/ReluRelu7/BiasAdd:output:0*
T0*'
_output_shapes
:���������@x
4/MatMul/ReadVariableOpReadVariableOp matmul_readvariableop_resource_1*
_output_shapes

:@@*
dtype0{
4/MatMulMatMulx/Relu:activations:04/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@v
4/BiasAdd/ReadVariableOpReadVariableOp!biasadd_readvariableop_resource_1*
_output_shapes
:@*
dtype0|
	4/BiasAddBiasAdd4/MatMul:product:0 4/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@T
4/ReluRelu4/BiasAdd:output:0*
T0*'
_output_shapes
:���������@x
1/MatMul/ReadVariableOpReadVariableOp matmul_readvariableop_resource_2*
_output_shapes

:@@*
dtype0{
1/MatMulMatMulx/Relu:activations:01/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@v
1/BiasAdd/ReadVariableOpReadVariableOp!biasadd_readvariableop_resource_2*
_output_shapes
:@*
dtype0|
	1/BiasAddBiasAdd1/MatMul:product:0 1/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@T
1/ReluRelu1/BiasAdd:output:0*
T0*'
_output_shapes
:���������@y
11/MatMul/ReadVariableOpReadVariableOp matmul_readvariableop_resource_3*
_output_shapes

:@@*
dtype0~
	11/MatMulMatMul10/Relu:activations:0 11/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@w
11/BiasAdd/ReadVariableOpReadVariableOp!biasadd_readvariableop_resource_3*
_output_shapes
:@*
dtype0

11/BiasAddBiasAdd11/MatMul:product:0!11/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@V
11/ReluRelu11/BiasAdd:output:0*
T0*'
_output_shapes
:���������@x
8/MatMul/ReadVariableOpReadVariableOp matmul_readvariableop_resource_4*
_output_shapes

:@@*
dtype0{
8/MatMulMatMul7/Relu:activations:08/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@v
8/BiasAdd/ReadVariableOpReadVariableOp!biasadd_readvariableop_resource_4*
_output_shapes
:@*
dtype0|
	8/BiasAddBiasAdd8/MatMul:product:0 8/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@T
8/ReluRelu8/BiasAdd:output:0*
T0*'
_output_shapes
:���������@x
5/MatMul/ReadVariableOpReadVariableOp matmul_readvariableop_resource_5*
_output_shapes

:@@*
dtype0{
5/MatMulMatMul4/Relu:activations:05/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@v
5/BiasAdd/ReadVariableOpReadVariableOp!biasadd_readvariableop_resource_5*
_output_shapes
:@*
dtype0|
	5/BiasAddBiasAdd5/MatMul:product:0 5/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@T
5/ReluRelu5/BiasAdd:output:0*
T0*'
_output_shapes
:���������@x
2/MatMul/ReadVariableOpReadVariableOp matmul_readvariableop_resource_6*
_output_shapes

:@@*
dtype0{
2/MatMulMatMul1/Relu:activations:02/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@v
2/BiasAdd/ReadVariableOpReadVariableOp!biasadd_readvariableop_resource_6*
_output_shapes
:@*
dtype0|
	2/BiasAddBiasAdd2/MatMul:product:0 2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@T
2/ReluRelu2/BiasAdd:output:0*
T0*'
_output_shapes
:���������@y
12/MatMul/ReadVariableOpReadVariableOp matmul_readvariableop_resource_7*
_output_shapes

:@ *
dtype0~
	12/MatMulMatMul11/Relu:activations:0 12/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� w
12/BiasAdd/ReadVariableOpReadVariableOp!biasadd_readvariableop_resource_7*
_output_shapes
: *
dtype0

12/BiasAddBiasAdd12/MatMul:product:0!12/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� V
12/ReluRelu12/BiasAdd:output:0*
T0*'
_output_shapes
:��������� x
9/MatMul/ReadVariableOpReadVariableOp matmul_readvariableop_resource_8*
_output_shapes

:@ *
dtype0{
9/MatMulMatMul8/Relu:activations:09/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� v
9/BiasAdd/ReadVariableOpReadVariableOp!biasadd_readvariableop_resource_8*
_output_shapes
: *
dtype0|
	9/BiasAddBiasAdd9/MatMul:product:0 9/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� T
9/ReluRelu9/BiasAdd:output:0*
T0*'
_output_shapes
:��������� x
6/MatMul/ReadVariableOpReadVariableOp matmul_readvariableop_resource_9*
_output_shapes

:@ *
dtype0{
6/MatMulMatMul5/Relu:activations:06/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� v
6/BiasAdd/ReadVariableOpReadVariableOp!biasadd_readvariableop_resource_9*
_output_shapes
: *
dtype0|
	6/BiasAddBiasAdd6/MatMul:product:0 6/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� T
6/ReluRelu6/BiasAdd:output:0*
T0*'
_output_shapes
:��������� y
3/MatMul/ReadVariableOpReadVariableOp!matmul_readvariableop_resource_10*
_output_shapes

:@ *
dtype0{
3/MatMulMatMul2/Relu:activations:03/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� w
3/BiasAdd/ReadVariableOpReadVariableOp"biasadd_readvariableop_resource_10*
_output_shapes
: *
dtype0|
	3/BiasAddBiasAdd3/MatMul:product:0 3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� T
3/ReluRelu3/BiasAdd:output:0*
T0*'
_output_shapes
:��������� x
c/MatMul/ReadVariableOpReadVariableOp c_matmul_readvariableop_resource*
_output_shapes

: *
dtype0|
c/MatMulMatMul12/Relu:activations:0c/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������v
c/BiasAdd/ReadVariableOpReadVariableOp!c_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0|
	c/BiasAddBiasAddc/MatMul:product:0 c/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������T
c/ReluReluc/BiasAdd:output:0*
T0*'
_output_shapes
:���������|
phi/MatMul/ReadVariableOpReadVariableOp"phi_matmul_readvariableop_resource*
_output_shapes

: *
dtype0

phi/MatMulMatMul9/Relu:activations:0!phi/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������z
phi/BiasAdd/ReadVariableOpReadVariableOp#phi_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
phi/BiasAddBiasAddphi/MatMul:product:0"phi/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������X
phi/ReluReluphi/BiasAdd:output:0*
T0*'
_output_shapes
:���������x
P/MatMul/ReadVariableOpReadVariableOp p_matmul_readvariableop_resource*
_output_shapes

: *
dtype0{
P/MatMulMatMul6/Relu:activations:0P/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������v
P/BiasAdd/ReadVariableOpReadVariableOp!p_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0|
	P/BiasAddBiasAddP/MatMul:product:0 P/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������T
P/ReluReluP/BiasAdd:output:0*
T0*'
_output_shapes
:���������x
T/MatMul/ReadVariableOpReadVariableOp t_matmul_readvariableop_resource*
_output_shapes

: *
dtype0{
T/MatMulMatMul3/Relu:activations:0T/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������v
T/BiasAdd/ReadVariableOpReadVariableOp!t_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0|
	T/BiasAddBiasAddT/MatMul:product:0 T/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������T
T/ReluReluT/BiasAdd:output:0*
T0*'
_output_shapes
:���������c
IdentityIdentityT/Relu:activations:0^NoOp*
T0*'
_output_shapes
:���������e

Identity_1IdentityP/Relu:activations:0^NoOp*
T0*'
_output_shapes
:���������g

Identity_2Identityphi/Relu:activations:0^NoOp*
T0*'
_output_shapes
:���������e

Identity_3Identityc/Relu:activations:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp^1/BiasAdd/ReadVariableOp^1/MatMul/ReadVariableOp^10/BiasAdd/ReadVariableOp^10/MatMul/ReadVariableOp^11/BiasAdd/ReadVariableOp^11/MatMul/ReadVariableOp^12/BiasAdd/ReadVariableOp^12/MatMul/ReadVariableOp^2/BiasAdd/ReadVariableOp^2/MatMul/ReadVariableOp^3/BiasAdd/ReadVariableOp^3/MatMul/ReadVariableOp^4/BiasAdd/ReadVariableOp^4/MatMul/ReadVariableOp^5/BiasAdd/ReadVariableOp^5/MatMul/ReadVariableOp^6/BiasAdd/ReadVariableOp^6/MatMul/ReadVariableOp^7/BiasAdd/ReadVariableOp^7/MatMul/ReadVariableOp^8/BiasAdd/ReadVariableOp^8/MatMul/ReadVariableOp^9/BiasAdd/ReadVariableOp^9/MatMul/ReadVariableOp^P/BiasAdd/ReadVariableOp^P/MatMul/ReadVariableOp^T/BiasAdd/ReadVariableOp^T/MatMul/ReadVariableOp^c/BiasAdd/ReadVariableOp^c/MatMul/ReadVariableOp^l2/BiasAdd/ReadVariableOp^l2/MatMul/ReadVariableOp^l3/BiasAdd/ReadVariableOp^l3/MatMul/ReadVariableOp^phi/BiasAdd/ReadVariableOp^phi/MatMul/ReadVariableOp^x/BiasAdd/ReadVariableOp^x/MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0"!

identity_3Identity_3:output:0*(
_construction_contextkEagerRuntime*�
_input_shapesu
s:���������::: : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : 24
1/BiasAdd/ReadVariableOp1/BiasAdd/ReadVariableOp22
1/MatMul/ReadVariableOp1/MatMul/ReadVariableOp26
10/BiasAdd/ReadVariableOp10/BiasAdd/ReadVariableOp24
10/MatMul/ReadVariableOp10/MatMul/ReadVariableOp26
11/BiasAdd/ReadVariableOp11/BiasAdd/ReadVariableOp24
11/MatMul/ReadVariableOp11/MatMul/ReadVariableOp26
12/BiasAdd/ReadVariableOp12/BiasAdd/ReadVariableOp24
12/MatMul/ReadVariableOp12/MatMul/ReadVariableOp24
2/BiasAdd/ReadVariableOp2/BiasAdd/ReadVariableOp22
2/MatMul/ReadVariableOp2/MatMul/ReadVariableOp24
3/BiasAdd/ReadVariableOp3/BiasAdd/ReadVariableOp22
3/MatMul/ReadVariableOp3/MatMul/ReadVariableOp24
4/BiasAdd/ReadVariableOp4/BiasAdd/ReadVariableOp22
4/MatMul/ReadVariableOp4/MatMul/ReadVariableOp24
5/BiasAdd/ReadVariableOp5/BiasAdd/ReadVariableOp22
5/MatMul/ReadVariableOp5/MatMul/ReadVariableOp24
6/BiasAdd/ReadVariableOp6/BiasAdd/ReadVariableOp22
6/MatMul/ReadVariableOp6/MatMul/ReadVariableOp24
7/BiasAdd/ReadVariableOp7/BiasAdd/ReadVariableOp22
7/MatMul/ReadVariableOp7/MatMul/ReadVariableOp24
8/BiasAdd/ReadVariableOp8/BiasAdd/ReadVariableOp22
8/MatMul/ReadVariableOp8/MatMul/ReadVariableOp24
9/BiasAdd/ReadVariableOp9/BiasAdd/ReadVariableOp22
9/MatMul/ReadVariableOp9/MatMul/ReadVariableOp24
P/BiasAdd/ReadVariableOpP/BiasAdd/ReadVariableOp22
P/MatMul/ReadVariableOpP/MatMul/ReadVariableOp24
T/BiasAdd/ReadVariableOpT/BiasAdd/ReadVariableOp22
T/MatMul/ReadVariableOpT/MatMul/ReadVariableOp24
c/BiasAdd/ReadVariableOpc/BiasAdd/ReadVariableOp22
c/MatMul/ReadVariableOpc/MatMul/ReadVariableOp26
l2/BiasAdd/ReadVariableOpl2/BiasAdd/ReadVariableOp24
l2/MatMul/ReadVariableOpl2/MatMul/ReadVariableOp26
l3/BiasAdd/ReadVariableOpl3/BiasAdd/ReadVariableOp24
l3/MatMul/ReadVariableOpl3/MatMul/ReadVariableOp28
phi/BiasAdd/ReadVariableOpphi/BiasAdd/ReadVariableOp26
phi/MatMul/ReadVariableOpphi/MatMul/ReadVariableOp24
x/BiasAdd/ReadVariableOpx/BiasAdd/ReadVariableOp22
x/MatMul/ReadVariableOpx/MatMul/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs:$ 

_output_shapes

::$ 

_output_shapes

:
�

�
>__inference_2_layer_call_and_return_conditional_losses_2217059

inputs0
matmul_readvariableop_resource:@@-
biasadd_readvariableop_resource:@
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@@*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������@a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:���������@w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�

�
>__inference_9_layer_call_and_return_conditional_losses_2218980

inputs0
matmul_readvariableop_resource:@ -
biasadd_readvariableop_resource: 
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@ *
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:��������� a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:��������� w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�

�
>__inference_7_layer_call_and_return_conditional_losses_2218820

inputs0
matmul_readvariableop_resource:@@-
biasadd_readvariableop_resource:@
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@@*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������@a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:���������@w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�

�
>__inference_8_layer_call_and_return_conditional_losses_2217025

inputs0
matmul_readvariableop_resource:@@-
biasadd_readvariableop_resource:@
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@@*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������@a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:���������@w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�
�
#__inference_T_layer_call_fn_2219009

inputs
unknown: 
	unknown_0:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_T_layer_call_and_return_conditional_losses_2217195o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:��������� : : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:��������� 
 
_user_specified_nameinputs
�

�
>__inference_x_layer_call_and_return_conditional_losses_2216923

inputs0
matmul_readvariableop_resource:@@-
biasadd_readvariableop_resource:@
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@@*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������@a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:���������@w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�

�
>__inference_1_layer_call_and_return_conditional_losses_2216991

inputs0
matmul_readvariableop_resource:@@-
biasadd_readvariableop_resource:@
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@@*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������@a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:���������@w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�
�
#__inference_9_layer_call_fn_2218969

inputs
unknown:@ 
	unknown_0: 
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:��������� *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_9_layer_call_and_return_conditional_losses_2217093o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:��������� `
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�

�
>__inference_6_layer_call_and_return_conditional_losses_2217110

inputs0
matmul_readvariableop_resource:@ -
biasadd_readvariableop_resource: 
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@ *
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:��������� a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:��������� w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
��
�
B__inference_model_layer_call_and_return_conditional_losses_2218700

inputs
normalization_sub_y
normalization_sqrt_x3
!l2_matmul_readvariableop_resource:@0
"l2_biasadd_readvariableop_resource:@3
!l3_matmul_readvariableop_resource:@@0
"l3_biasadd_readvariableop_resource:@2
 x_matmul_readvariableop_resource:@@/
!x_biasadd_readvariableop_resource:@0
matmul_readvariableop_resource:@@-
biasadd_readvariableop_resource:@2
 matmul_readvariableop_resource_0:@@/
!biasadd_readvariableop_resource_0:@2
 matmul_readvariableop_resource_1:@@/
!biasadd_readvariableop_resource_1:@2
 matmul_readvariableop_resource_2:@@/
!biasadd_readvariableop_resource_2:@2
 matmul_readvariableop_resource_3:@@/
!biasadd_readvariableop_resource_3:@2
 matmul_readvariableop_resource_4:@@/
!biasadd_readvariableop_resource_4:@2
 matmul_readvariableop_resource_5:@@/
!biasadd_readvariableop_resource_5:@2
 matmul_readvariableop_resource_6:@@/
!biasadd_readvariableop_resource_6:@2
 matmul_readvariableop_resource_7:@ /
!biasadd_readvariableop_resource_7: 2
 matmul_readvariableop_resource_8:@ /
!biasadd_readvariableop_resource_8: 2
 matmul_readvariableop_resource_9:@ /
!biasadd_readvariableop_resource_9: 3
!matmul_readvariableop_resource_10:@ 0
"biasadd_readvariableop_resource_10: 2
 c_matmul_readvariableop_resource: /
!c_biasadd_readvariableop_resource:4
"phi_matmul_readvariableop_resource: 1
#phi_biasadd_readvariableop_resource:2
 p_matmul_readvariableop_resource: /
!p_biasadd_readvariableop_resource:2
 t_matmul_readvariableop_resource: /
!t_biasadd_readvariableop_resource:
identity

identity_1

identity_2

identity_3��1/BiasAdd/ReadVariableOp�1/MatMul/ReadVariableOp�10/BiasAdd/ReadVariableOp�10/MatMul/ReadVariableOp�11/BiasAdd/ReadVariableOp�11/MatMul/ReadVariableOp�12/BiasAdd/ReadVariableOp�12/MatMul/ReadVariableOp�2/BiasAdd/ReadVariableOp�2/MatMul/ReadVariableOp�3/BiasAdd/ReadVariableOp�3/MatMul/ReadVariableOp�4/BiasAdd/ReadVariableOp�4/MatMul/ReadVariableOp�5/BiasAdd/ReadVariableOp�5/MatMul/ReadVariableOp�6/BiasAdd/ReadVariableOp�6/MatMul/ReadVariableOp�7/BiasAdd/ReadVariableOp�7/MatMul/ReadVariableOp�8/BiasAdd/ReadVariableOp�8/MatMul/ReadVariableOp�9/BiasAdd/ReadVariableOp�9/MatMul/ReadVariableOp�P/BiasAdd/ReadVariableOp�P/MatMul/ReadVariableOp�T/BiasAdd/ReadVariableOp�T/MatMul/ReadVariableOp�c/BiasAdd/ReadVariableOp�c/MatMul/ReadVariableOp�l2/BiasAdd/ReadVariableOp�l2/MatMul/ReadVariableOp�l3/BiasAdd/ReadVariableOp�l3/MatMul/ReadVariableOp�phi/BiasAdd/ReadVariableOp�phi/MatMul/ReadVariableOp�x/BiasAdd/ReadVariableOp�x/MatMul/ReadVariableOpg
normalization/subSubinputsnormalization_sub_y*
T0*'
_output_shapes
:���������Y
normalization/SqrtSqrtnormalization_sqrt_x*
T0*
_output_shapes

:\
normalization/Maximum/yConst*
_output_shapes
: *
dtype0*
valueB
 *���3�
normalization/MaximumMaximumnormalization/Sqrt:y:0 normalization/Maximum/y:output:0*
T0*
_output_shapes

:�
normalization/truedivRealDivnormalization/sub:z:0normalization/Maximum:z:0*
T0*'
_output_shapes
:���������z
l2/MatMul/ReadVariableOpReadVariableOp!l2_matmul_readvariableop_resource*
_output_shapes

:@*
dtype0�
	l2/MatMulMatMulnormalization/truediv:z:0 l2/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@x
l2/BiasAdd/ReadVariableOpReadVariableOp"l2_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype0

l2/BiasAddBiasAddl2/MatMul:product:0!l2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@V
l2/ReluRelul2/BiasAdd:output:0*
T0*'
_output_shapes
:���������@z
l3/MatMul/ReadVariableOpReadVariableOp!l3_matmul_readvariableop_resource*
_output_shapes

:@@*
dtype0~
	l3/MatMulMatMull2/Relu:activations:0 l3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@x
l3/BiasAdd/ReadVariableOpReadVariableOp"l3_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype0

l3/BiasAddBiasAddl3/MatMul:product:0!l3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@V
l3/ReluRelul3/BiasAdd:output:0*
T0*'
_output_shapes
:���������@x
x/MatMul/ReadVariableOpReadVariableOp x_matmul_readvariableop_resource*
_output_shapes

:@@*
dtype0|
x/MatMulMatMull3/Relu:activations:0x/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@v
x/BiasAdd/ReadVariableOpReadVariableOp!x_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype0|
	x/BiasAddBiasAddx/MatMul:product:0 x/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@T
x/ReluRelux/BiasAdd:output:0*
T0*'
_output_shapes
:���������@w
10/MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@@*
dtype0}
	10/MatMulMatMulx/Relu:activations:0 10/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@u
10/BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype0

10/BiasAddBiasAdd10/MatMul:product:0!10/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@V
10/ReluRelu10/BiasAdd:output:0*
T0*'
_output_shapes
:���������@x
7/MatMul/ReadVariableOpReadVariableOp matmul_readvariableop_resource_0*
_output_shapes

:@@*
dtype0{
7/MatMulMatMulx/Relu:activations:07/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@v
7/BiasAdd/ReadVariableOpReadVariableOp!biasadd_readvariableop_resource_0*
_output_shapes
:@*
dtype0|
	7/BiasAddBiasAdd7/MatMul:product:0 7/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@T
7/ReluRelu7/BiasAdd:output:0*
T0*'
_output_shapes
:���������@x
4/MatMul/ReadVariableOpReadVariableOp matmul_readvariableop_resource_1*
_output_shapes

:@@*
dtype0{
4/MatMulMatMulx/Relu:activations:04/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@v
4/BiasAdd/ReadVariableOpReadVariableOp!biasadd_readvariableop_resource_1*
_output_shapes
:@*
dtype0|
	4/BiasAddBiasAdd4/MatMul:product:0 4/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@T
4/ReluRelu4/BiasAdd:output:0*
T0*'
_output_shapes
:���������@x
1/MatMul/ReadVariableOpReadVariableOp matmul_readvariableop_resource_2*
_output_shapes

:@@*
dtype0{
1/MatMulMatMulx/Relu:activations:01/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@v
1/BiasAdd/ReadVariableOpReadVariableOp!biasadd_readvariableop_resource_2*
_output_shapes
:@*
dtype0|
	1/BiasAddBiasAdd1/MatMul:product:0 1/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@T
1/ReluRelu1/BiasAdd:output:0*
T0*'
_output_shapes
:���������@y
11/MatMul/ReadVariableOpReadVariableOp matmul_readvariableop_resource_3*
_output_shapes

:@@*
dtype0~
	11/MatMulMatMul10/Relu:activations:0 11/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@w
11/BiasAdd/ReadVariableOpReadVariableOp!biasadd_readvariableop_resource_3*
_output_shapes
:@*
dtype0

11/BiasAddBiasAdd11/MatMul:product:0!11/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@V
11/ReluRelu11/BiasAdd:output:0*
T0*'
_output_shapes
:���������@x
8/MatMul/ReadVariableOpReadVariableOp matmul_readvariableop_resource_4*
_output_shapes

:@@*
dtype0{
8/MatMulMatMul7/Relu:activations:08/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@v
8/BiasAdd/ReadVariableOpReadVariableOp!biasadd_readvariableop_resource_4*
_output_shapes
:@*
dtype0|
	8/BiasAddBiasAdd8/MatMul:product:0 8/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@T
8/ReluRelu8/BiasAdd:output:0*
T0*'
_output_shapes
:���������@x
5/MatMul/ReadVariableOpReadVariableOp matmul_readvariableop_resource_5*
_output_shapes

:@@*
dtype0{
5/MatMulMatMul4/Relu:activations:05/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@v
5/BiasAdd/ReadVariableOpReadVariableOp!biasadd_readvariableop_resource_5*
_output_shapes
:@*
dtype0|
	5/BiasAddBiasAdd5/MatMul:product:0 5/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@T
5/ReluRelu5/BiasAdd:output:0*
T0*'
_output_shapes
:���������@x
2/MatMul/ReadVariableOpReadVariableOp matmul_readvariableop_resource_6*
_output_shapes

:@@*
dtype0{
2/MatMulMatMul1/Relu:activations:02/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@v
2/BiasAdd/ReadVariableOpReadVariableOp!biasadd_readvariableop_resource_6*
_output_shapes
:@*
dtype0|
	2/BiasAddBiasAdd2/MatMul:product:0 2/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@T
2/ReluRelu2/BiasAdd:output:0*
T0*'
_output_shapes
:���������@y
12/MatMul/ReadVariableOpReadVariableOp matmul_readvariableop_resource_7*
_output_shapes

:@ *
dtype0~
	12/MatMulMatMul11/Relu:activations:0 12/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� w
12/BiasAdd/ReadVariableOpReadVariableOp!biasadd_readvariableop_resource_7*
_output_shapes
: *
dtype0

12/BiasAddBiasAdd12/MatMul:product:0!12/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� V
12/ReluRelu12/BiasAdd:output:0*
T0*'
_output_shapes
:��������� x
9/MatMul/ReadVariableOpReadVariableOp matmul_readvariableop_resource_8*
_output_shapes

:@ *
dtype0{
9/MatMulMatMul8/Relu:activations:09/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� v
9/BiasAdd/ReadVariableOpReadVariableOp!biasadd_readvariableop_resource_8*
_output_shapes
: *
dtype0|
	9/BiasAddBiasAdd9/MatMul:product:0 9/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� T
9/ReluRelu9/BiasAdd:output:0*
T0*'
_output_shapes
:��������� x
6/MatMul/ReadVariableOpReadVariableOp matmul_readvariableop_resource_9*
_output_shapes

:@ *
dtype0{
6/MatMulMatMul5/Relu:activations:06/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� v
6/BiasAdd/ReadVariableOpReadVariableOp!biasadd_readvariableop_resource_9*
_output_shapes
: *
dtype0|
	6/BiasAddBiasAdd6/MatMul:product:0 6/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� T
6/ReluRelu6/BiasAdd:output:0*
T0*'
_output_shapes
:��������� y
3/MatMul/ReadVariableOpReadVariableOp!matmul_readvariableop_resource_10*
_output_shapes

:@ *
dtype0{
3/MatMulMatMul2/Relu:activations:03/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� w
3/BiasAdd/ReadVariableOpReadVariableOp"biasadd_readvariableop_resource_10*
_output_shapes
: *
dtype0|
	3/BiasAddBiasAdd3/MatMul:product:0 3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� T
3/ReluRelu3/BiasAdd:output:0*
T0*'
_output_shapes
:��������� x
c/MatMul/ReadVariableOpReadVariableOp c_matmul_readvariableop_resource*
_output_shapes

: *
dtype0|
c/MatMulMatMul12/Relu:activations:0c/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������v
c/BiasAdd/ReadVariableOpReadVariableOp!c_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0|
	c/BiasAddBiasAddc/MatMul:product:0 c/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������T
c/ReluReluc/BiasAdd:output:0*
T0*'
_output_shapes
:���������|
phi/MatMul/ReadVariableOpReadVariableOp"phi_matmul_readvariableop_resource*
_output_shapes

: *
dtype0

phi/MatMulMatMul9/Relu:activations:0!phi/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������z
phi/BiasAdd/ReadVariableOpReadVariableOp#phi_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0�
phi/BiasAddBiasAddphi/MatMul:product:0"phi/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������X
phi/ReluReluphi/BiasAdd:output:0*
T0*'
_output_shapes
:���������x
P/MatMul/ReadVariableOpReadVariableOp p_matmul_readvariableop_resource*
_output_shapes

: *
dtype0{
P/MatMulMatMul6/Relu:activations:0P/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������v
P/BiasAdd/ReadVariableOpReadVariableOp!p_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0|
	P/BiasAddBiasAddP/MatMul:product:0 P/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������T
P/ReluReluP/BiasAdd:output:0*
T0*'
_output_shapes
:���������x
T/MatMul/ReadVariableOpReadVariableOp t_matmul_readvariableop_resource*
_output_shapes

: *
dtype0{
T/MatMulMatMul3/Relu:activations:0T/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������v
T/BiasAdd/ReadVariableOpReadVariableOp!t_biasadd_readvariableop_resource*
_output_shapes
:*
dtype0|
	T/BiasAddBiasAddT/MatMul:product:0 T/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������T
T/ReluReluT/BiasAdd:output:0*
T0*'
_output_shapes
:���������c
IdentityIdentityT/Relu:activations:0^NoOp*
T0*'
_output_shapes
:���������e

Identity_1IdentityP/Relu:activations:0^NoOp*
T0*'
_output_shapes
:���������g

Identity_2Identityphi/Relu:activations:0^NoOp*
T0*'
_output_shapes
:���������e

Identity_3Identityc/Relu:activations:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp^1/BiasAdd/ReadVariableOp^1/MatMul/ReadVariableOp^10/BiasAdd/ReadVariableOp^10/MatMul/ReadVariableOp^11/BiasAdd/ReadVariableOp^11/MatMul/ReadVariableOp^12/BiasAdd/ReadVariableOp^12/MatMul/ReadVariableOp^2/BiasAdd/ReadVariableOp^2/MatMul/ReadVariableOp^3/BiasAdd/ReadVariableOp^3/MatMul/ReadVariableOp^4/BiasAdd/ReadVariableOp^4/MatMul/ReadVariableOp^5/BiasAdd/ReadVariableOp^5/MatMul/ReadVariableOp^6/BiasAdd/ReadVariableOp^6/MatMul/ReadVariableOp^7/BiasAdd/ReadVariableOp^7/MatMul/ReadVariableOp^8/BiasAdd/ReadVariableOp^8/MatMul/ReadVariableOp^9/BiasAdd/ReadVariableOp^9/MatMul/ReadVariableOp^P/BiasAdd/ReadVariableOp^P/MatMul/ReadVariableOp^T/BiasAdd/ReadVariableOp^T/MatMul/ReadVariableOp^c/BiasAdd/ReadVariableOp^c/MatMul/ReadVariableOp^l2/BiasAdd/ReadVariableOp^l2/MatMul/ReadVariableOp^l3/BiasAdd/ReadVariableOp^l3/MatMul/ReadVariableOp^phi/BiasAdd/ReadVariableOp^phi/MatMul/ReadVariableOp^x/BiasAdd/ReadVariableOp^x/MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0"!

identity_3Identity_3:output:0*(
_construction_contextkEagerRuntime*�
_input_shapesu
s:���������::: : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : 24
1/BiasAdd/ReadVariableOp1/BiasAdd/ReadVariableOp22
1/MatMul/ReadVariableOp1/MatMul/ReadVariableOp26
10/BiasAdd/ReadVariableOp10/BiasAdd/ReadVariableOp24
10/MatMul/ReadVariableOp10/MatMul/ReadVariableOp26
11/BiasAdd/ReadVariableOp11/BiasAdd/ReadVariableOp24
11/MatMul/ReadVariableOp11/MatMul/ReadVariableOp26
12/BiasAdd/ReadVariableOp12/BiasAdd/ReadVariableOp24
12/MatMul/ReadVariableOp12/MatMul/ReadVariableOp24
2/BiasAdd/ReadVariableOp2/BiasAdd/ReadVariableOp22
2/MatMul/ReadVariableOp2/MatMul/ReadVariableOp24
3/BiasAdd/ReadVariableOp3/BiasAdd/ReadVariableOp22
3/MatMul/ReadVariableOp3/MatMul/ReadVariableOp24
4/BiasAdd/ReadVariableOp4/BiasAdd/ReadVariableOp22
4/MatMul/ReadVariableOp4/MatMul/ReadVariableOp24
5/BiasAdd/ReadVariableOp5/BiasAdd/ReadVariableOp22
5/MatMul/ReadVariableOp5/MatMul/ReadVariableOp24
6/BiasAdd/ReadVariableOp6/BiasAdd/ReadVariableOp22
6/MatMul/ReadVariableOp6/MatMul/ReadVariableOp24
7/BiasAdd/ReadVariableOp7/BiasAdd/ReadVariableOp22
7/MatMul/ReadVariableOp7/MatMul/ReadVariableOp24
8/BiasAdd/ReadVariableOp8/BiasAdd/ReadVariableOp22
8/MatMul/ReadVariableOp8/MatMul/ReadVariableOp24
9/BiasAdd/ReadVariableOp9/BiasAdd/ReadVariableOp22
9/MatMul/ReadVariableOp9/MatMul/ReadVariableOp24
P/BiasAdd/ReadVariableOpP/BiasAdd/ReadVariableOp22
P/MatMul/ReadVariableOpP/MatMul/ReadVariableOp24
T/BiasAdd/ReadVariableOpT/BiasAdd/ReadVariableOp22
T/MatMul/ReadVariableOpT/MatMul/ReadVariableOp24
c/BiasAdd/ReadVariableOpc/BiasAdd/ReadVariableOp22
c/MatMul/ReadVariableOpc/MatMul/ReadVariableOp26
l2/BiasAdd/ReadVariableOpl2/BiasAdd/ReadVariableOp24
l2/MatMul/ReadVariableOpl2/MatMul/ReadVariableOp26
l3/BiasAdd/ReadVariableOpl3/BiasAdd/ReadVariableOp24
l3/MatMul/ReadVariableOpl3/MatMul/ReadVariableOp28
phi/BiasAdd/ReadVariableOpphi/BiasAdd/ReadVariableOp26
phi/MatMul/ReadVariableOpphi/MatMul/ReadVariableOp24
x/BiasAdd/ReadVariableOpx/BiasAdd/ReadVariableOp22
x/MatMul/ReadVariableOpx/MatMul/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs:$ 

_output_shapes

::$ 

_output_shapes

:
�
�H
#__inference__traced_restore_2219892
file_prefix#
assignvariableop_mean:)
assignvariableop_1_variance:$
assignvariableop_2_count_5:	 .
assignvariableop_3_l2_kernel:@(
assignvariableop_4_l2_bias:@.
assignvariableop_5_l3_kernel:@@(
assignvariableop_6_l3_bias:@-
assignvariableop_7_x_kernel:@@'
assignvariableop_8_x_bias:@-
assignvariableop_9_1_kernel:@@(
assignvariableop_10_1_bias:@.
assignvariableop_11_4_kernel:@@(
assignvariableop_12_4_bias:@.
assignvariableop_13_7_kernel:@@(
assignvariableop_14_7_bias:@/
assignvariableop_15_10_kernel:@@)
assignvariableop_16_10_bias:@.
assignvariableop_17_2_kernel:@@(
assignvariableop_18_2_bias:@.
assignvariableop_19_5_kernel:@@(
assignvariableop_20_5_bias:@.
assignvariableop_21_8_kernel:@@(
assignvariableop_22_8_bias:@/
assignvariableop_23_11_kernel:@@)
assignvariableop_24_11_bias:@.
assignvariableop_25_3_kernel:@ (
assignvariableop_26_3_bias: .
assignvariableop_27_6_kernel:@ (
assignvariableop_28_6_bias: .
assignvariableop_29_9_kernel:@ (
assignvariableop_30_9_bias: /
assignvariableop_31_12_kernel:@ )
assignvariableop_32_12_bias: .
assignvariableop_33_t_kernel: (
assignvariableop_34_t_bias:.
assignvariableop_35_p_kernel: (
assignvariableop_36_p_bias:0
assignvariableop_37_phi_kernel: *
assignvariableop_38_phi_bias:.
assignvariableop_39_c_kernel: (
assignvariableop_40_c_bias:'
assignvariableop_41_iteration:	 +
!assignvariableop_42_learning_rate: 6
$assignvariableop_43_adam_m_l2_kernel:@6
$assignvariableop_44_adam_v_l2_kernel:@0
"assignvariableop_45_adam_m_l2_bias:@0
"assignvariableop_46_adam_v_l2_bias:@6
$assignvariableop_47_adam_m_l3_kernel:@@6
$assignvariableop_48_adam_v_l3_kernel:@@0
"assignvariableop_49_adam_m_l3_bias:@0
"assignvariableop_50_adam_v_l3_bias:@5
#assignvariableop_51_adam_m_x_kernel:@@5
#assignvariableop_52_adam_v_x_kernel:@@/
!assignvariableop_53_adam_m_x_bias:@/
!assignvariableop_54_adam_v_x_bias:@5
#assignvariableop_55_adam_m_1_kernel:@@5
#assignvariableop_56_adam_v_1_kernel:@@/
!assignvariableop_57_adam_m_1_bias:@/
!assignvariableop_58_adam_v_1_bias:@5
#assignvariableop_59_adam_m_4_kernel:@@5
#assignvariableop_60_adam_v_4_kernel:@@/
!assignvariableop_61_adam_m_4_bias:@/
!assignvariableop_62_adam_v_4_bias:@5
#assignvariableop_63_adam_m_7_kernel:@@5
#assignvariableop_64_adam_v_7_kernel:@@/
!assignvariableop_65_adam_m_7_bias:@/
!assignvariableop_66_adam_v_7_bias:@6
$assignvariableop_67_adam_m_10_kernel:@@6
$assignvariableop_68_adam_v_10_kernel:@@0
"assignvariableop_69_adam_m_10_bias:@0
"assignvariableop_70_adam_v_10_bias:@5
#assignvariableop_71_adam_m_2_kernel:@@5
#assignvariableop_72_adam_v_2_kernel:@@/
!assignvariableop_73_adam_m_2_bias:@/
!assignvariableop_74_adam_v_2_bias:@5
#assignvariableop_75_adam_m_5_kernel:@@5
#assignvariableop_76_adam_v_5_kernel:@@/
!assignvariableop_77_adam_m_5_bias:@/
!assignvariableop_78_adam_v_5_bias:@5
#assignvariableop_79_adam_m_8_kernel:@@5
#assignvariableop_80_adam_v_8_kernel:@@/
!assignvariableop_81_adam_m_8_bias:@/
!assignvariableop_82_adam_v_8_bias:@6
$assignvariableop_83_adam_m_11_kernel:@@6
$assignvariableop_84_adam_v_11_kernel:@@0
"assignvariableop_85_adam_m_11_bias:@0
"assignvariableop_86_adam_v_11_bias:@5
#assignvariableop_87_adam_m_3_kernel:@ 5
#assignvariableop_88_adam_v_3_kernel:@ /
!assignvariableop_89_adam_m_3_bias: /
!assignvariableop_90_adam_v_3_bias: 5
#assignvariableop_91_adam_m_6_kernel:@ 5
#assignvariableop_92_adam_v_6_kernel:@ /
!assignvariableop_93_adam_m_6_bias: /
!assignvariableop_94_adam_v_6_bias: 5
#assignvariableop_95_adam_m_9_kernel:@ 5
#assignvariableop_96_adam_v_9_kernel:@ /
!assignvariableop_97_adam_m_9_bias: /
!assignvariableop_98_adam_v_9_bias: 6
$assignvariableop_99_adam_m_12_kernel:@ 7
%assignvariableop_100_adam_v_12_kernel:@ 1
#assignvariableop_101_adam_m_12_bias: 1
#assignvariableop_102_adam_v_12_bias: 6
$assignvariableop_103_adam_m_t_kernel: 6
$assignvariableop_104_adam_v_t_kernel: 0
"assignvariableop_105_adam_m_t_bias:0
"assignvariableop_106_adam_v_t_bias:6
$assignvariableop_107_adam_m_p_kernel: 6
$assignvariableop_108_adam_v_p_kernel: 0
"assignvariableop_109_adam_m_p_bias:0
"assignvariableop_110_adam_v_p_bias:8
&assignvariableop_111_adam_m_phi_kernel: 8
&assignvariableop_112_adam_v_phi_kernel: 2
$assignvariableop_113_adam_m_phi_bias:2
$assignvariableop_114_adam_v_phi_bias:6
$assignvariableop_115_adam_m_c_kernel: 6
$assignvariableop_116_adam_v_c_kernel: 0
"assignvariableop_117_adam_m_c_bias:0
"assignvariableop_118_adam_v_c_bias:&
assignvariableop_119_total_4: &
assignvariableop_120_count_4: &
assignvariableop_121_total_3: &
assignvariableop_122_count_3: &
assignvariableop_123_total_2: &
assignvariableop_124_count_2: &
assignvariableop_125_total_1: &
assignvariableop_126_count_1: $
assignvariableop_127_total: $
assignvariableop_128_count: 
identity_130��AssignVariableOp�AssignVariableOp_1�AssignVariableOp_10�AssignVariableOp_100�AssignVariableOp_101�AssignVariableOp_102�AssignVariableOp_103�AssignVariableOp_104�AssignVariableOp_105�AssignVariableOp_106�AssignVariableOp_107�AssignVariableOp_108�AssignVariableOp_109�AssignVariableOp_11�AssignVariableOp_110�AssignVariableOp_111�AssignVariableOp_112�AssignVariableOp_113�AssignVariableOp_114�AssignVariableOp_115�AssignVariableOp_116�AssignVariableOp_117�AssignVariableOp_118�AssignVariableOp_119�AssignVariableOp_12�AssignVariableOp_120�AssignVariableOp_121�AssignVariableOp_122�AssignVariableOp_123�AssignVariableOp_124�AssignVariableOp_125�AssignVariableOp_126�AssignVariableOp_127�AssignVariableOp_128�AssignVariableOp_13�AssignVariableOp_14�AssignVariableOp_15�AssignVariableOp_16�AssignVariableOp_17�AssignVariableOp_18�AssignVariableOp_19�AssignVariableOp_2�AssignVariableOp_20�AssignVariableOp_21�AssignVariableOp_22�AssignVariableOp_23�AssignVariableOp_24�AssignVariableOp_25�AssignVariableOp_26�AssignVariableOp_27�AssignVariableOp_28�AssignVariableOp_29�AssignVariableOp_3�AssignVariableOp_30�AssignVariableOp_31�AssignVariableOp_32�AssignVariableOp_33�AssignVariableOp_34�AssignVariableOp_35�AssignVariableOp_36�AssignVariableOp_37�AssignVariableOp_38�AssignVariableOp_39�AssignVariableOp_4�AssignVariableOp_40�AssignVariableOp_41�AssignVariableOp_42�AssignVariableOp_43�AssignVariableOp_44�AssignVariableOp_45�AssignVariableOp_46�AssignVariableOp_47�AssignVariableOp_48�AssignVariableOp_49�AssignVariableOp_5�AssignVariableOp_50�AssignVariableOp_51�AssignVariableOp_52�AssignVariableOp_53�AssignVariableOp_54�AssignVariableOp_55�AssignVariableOp_56�AssignVariableOp_57�AssignVariableOp_58�AssignVariableOp_59�AssignVariableOp_6�AssignVariableOp_60�AssignVariableOp_61�AssignVariableOp_62�AssignVariableOp_63�AssignVariableOp_64�AssignVariableOp_65�AssignVariableOp_66�AssignVariableOp_67�AssignVariableOp_68�AssignVariableOp_69�AssignVariableOp_7�AssignVariableOp_70�AssignVariableOp_71�AssignVariableOp_72�AssignVariableOp_73�AssignVariableOp_74�AssignVariableOp_75�AssignVariableOp_76�AssignVariableOp_77�AssignVariableOp_78�AssignVariableOp_79�AssignVariableOp_8�AssignVariableOp_80�AssignVariableOp_81�AssignVariableOp_82�AssignVariableOp_83�AssignVariableOp_84�AssignVariableOp_85�AssignVariableOp_86�AssignVariableOp_87�AssignVariableOp_88�AssignVariableOp_89�AssignVariableOp_9�AssignVariableOp_90�AssignVariableOp_91�AssignVariableOp_92�AssignVariableOp_93�AssignVariableOp_94�AssignVariableOp_95�AssignVariableOp_96�AssignVariableOp_97�AssignVariableOp_98�AssignVariableOp_99�6
RestoreV2/tensor_namesConst"/device:CPU:0*
_output_shapes	
:�*
dtype0*�6
value�5B�5�B4layer_with_weights-0/mean/.ATTRIBUTES/VARIABLE_VALUEB8layer_with_weights-0/variance/.ATTRIBUTES/VARIABLE_VALUEB5layer_with_weights-0/count/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-3/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-3/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-4/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-4/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-5/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-5/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-6/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-6/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-7/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-7/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-8/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-8/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-9/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-9/bias/.ATTRIBUTES/VARIABLE_VALUEB7layer_with_weights-10/kernel/.ATTRIBUTES/VARIABLE_VALUEB5layer_with_weights-10/bias/.ATTRIBUTES/VARIABLE_VALUEB7layer_with_weights-11/kernel/.ATTRIBUTES/VARIABLE_VALUEB5layer_with_weights-11/bias/.ATTRIBUTES/VARIABLE_VALUEB7layer_with_weights-12/kernel/.ATTRIBUTES/VARIABLE_VALUEB5layer_with_weights-12/bias/.ATTRIBUTES/VARIABLE_VALUEB7layer_with_weights-13/kernel/.ATTRIBUTES/VARIABLE_VALUEB5layer_with_weights-13/bias/.ATTRIBUTES/VARIABLE_VALUEB7layer_with_weights-14/kernel/.ATTRIBUTES/VARIABLE_VALUEB5layer_with_weights-14/bias/.ATTRIBUTES/VARIABLE_VALUEB7layer_with_weights-15/kernel/.ATTRIBUTES/VARIABLE_VALUEB5layer_with_weights-15/bias/.ATTRIBUTES/VARIABLE_VALUEB7layer_with_weights-16/kernel/.ATTRIBUTES/VARIABLE_VALUEB5layer_with_weights-16/bias/.ATTRIBUTES/VARIABLE_VALUEB7layer_with_weights-17/kernel/.ATTRIBUTES/VARIABLE_VALUEB5layer_with_weights-17/bias/.ATTRIBUTES/VARIABLE_VALUEB7layer_with_weights-18/kernel/.ATTRIBUTES/VARIABLE_VALUEB5layer_with_weights-18/bias/.ATTRIBUTES/VARIABLE_VALUEB7layer_with_weights-19/kernel/.ATTRIBUTES/VARIABLE_VALUEB5layer_with_weights-19/bias/.ATTRIBUTES/VARIABLE_VALUEB0optimizer/_iterations/.ATTRIBUTES/VARIABLE_VALUEB3optimizer/_learning_rate/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/1/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/2/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/3/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/4/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/5/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/6/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/7/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/8/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/9/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/10/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/11/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/12/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/13/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/14/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/15/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/16/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/17/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/18/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/19/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/20/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/21/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/22/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/23/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/24/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/25/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/26/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/27/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/28/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/29/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/30/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/31/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/32/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/33/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/34/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/35/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/36/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/37/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/38/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/39/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/40/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/41/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/42/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/43/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/44/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/45/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/46/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/47/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/48/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/49/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/50/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/51/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/52/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/53/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/54/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/55/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/56/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/57/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/58/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/59/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/60/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/61/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/62/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/63/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/64/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/65/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/66/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/67/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/68/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/69/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/70/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/71/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/72/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/73/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/74/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/75/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/76/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/2/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/2/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/3/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/3/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/4/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/4/count/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH�
RestoreV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes	
:�*
dtype0*�
value�B��B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B �
	RestoreV2	RestoreV2file_prefixRestoreV2/tensor_names:output:0#RestoreV2/shape_and_slices:output:0"/device:CPU:0*�
_output_shapes�
�::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*�
dtypes�
�2�		[
IdentityIdentityRestoreV2:tensors:0"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOpAssignVariableOpassignvariableop_meanIdentity:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_1IdentityRestoreV2:tensors:1"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_1AssignVariableOpassignvariableop_1_varianceIdentity_1:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_2IdentityRestoreV2:tensors:2"/device:CPU:0*
T0	*
_output_shapes
:�
AssignVariableOp_2AssignVariableOpassignvariableop_2_count_5Identity_2:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0	]

Identity_3IdentityRestoreV2:tensors:3"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_3AssignVariableOpassignvariableop_3_l2_kernelIdentity_3:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_4IdentityRestoreV2:tensors:4"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_4AssignVariableOpassignvariableop_4_l2_biasIdentity_4:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_5IdentityRestoreV2:tensors:5"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_5AssignVariableOpassignvariableop_5_l3_kernelIdentity_5:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_6IdentityRestoreV2:tensors:6"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_6AssignVariableOpassignvariableop_6_l3_biasIdentity_6:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_7IdentityRestoreV2:tensors:7"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_7AssignVariableOpassignvariableop_7_x_kernelIdentity_7:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_8IdentityRestoreV2:tensors:8"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_8AssignVariableOpassignvariableop_8_x_biasIdentity_8:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0]

Identity_9IdentityRestoreV2:tensors:9"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_9AssignVariableOpassignvariableop_9_1_kernelIdentity_9:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_10IdentityRestoreV2:tensors:10"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_10AssignVariableOpassignvariableop_10_1_biasIdentity_10:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_11IdentityRestoreV2:tensors:11"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_11AssignVariableOpassignvariableop_11_4_kernelIdentity_11:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_12IdentityRestoreV2:tensors:12"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_12AssignVariableOpassignvariableop_12_4_biasIdentity_12:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_13IdentityRestoreV2:tensors:13"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_13AssignVariableOpassignvariableop_13_7_kernelIdentity_13:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_14IdentityRestoreV2:tensors:14"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_14AssignVariableOpassignvariableop_14_7_biasIdentity_14:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_15IdentityRestoreV2:tensors:15"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_15AssignVariableOpassignvariableop_15_10_kernelIdentity_15:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_16IdentityRestoreV2:tensors:16"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_16AssignVariableOpassignvariableop_16_10_biasIdentity_16:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_17IdentityRestoreV2:tensors:17"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_17AssignVariableOpassignvariableop_17_2_kernelIdentity_17:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_18IdentityRestoreV2:tensors:18"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_18AssignVariableOpassignvariableop_18_2_biasIdentity_18:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_19IdentityRestoreV2:tensors:19"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_19AssignVariableOpassignvariableop_19_5_kernelIdentity_19:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_20IdentityRestoreV2:tensors:20"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_20AssignVariableOpassignvariableop_20_5_biasIdentity_20:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_21IdentityRestoreV2:tensors:21"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_21AssignVariableOpassignvariableop_21_8_kernelIdentity_21:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_22IdentityRestoreV2:tensors:22"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_22AssignVariableOpassignvariableop_22_8_biasIdentity_22:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_23IdentityRestoreV2:tensors:23"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_23AssignVariableOpassignvariableop_23_11_kernelIdentity_23:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_24IdentityRestoreV2:tensors:24"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_24AssignVariableOpassignvariableop_24_11_biasIdentity_24:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_25IdentityRestoreV2:tensors:25"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_25AssignVariableOpassignvariableop_25_3_kernelIdentity_25:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_26IdentityRestoreV2:tensors:26"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_26AssignVariableOpassignvariableop_26_3_biasIdentity_26:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_27IdentityRestoreV2:tensors:27"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_27AssignVariableOpassignvariableop_27_6_kernelIdentity_27:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_28IdentityRestoreV2:tensors:28"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_28AssignVariableOpassignvariableop_28_6_biasIdentity_28:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_29IdentityRestoreV2:tensors:29"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_29AssignVariableOpassignvariableop_29_9_kernelIdentity_29:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_30IdentityRestoreV2:tensors:30"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_30AssignVariableOpassignvariableop_30_9_biasIdentity_30:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_31IdentityRestoreV2:tensors:31"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_31AssignVariableOpassignvariableop_31_12_kernelIdentity_31:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_32IdentityRestoreV2:tensors:32"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_32AssignVariableOpassignvariableop_32_12_biasIdentity_32:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_33IdentityRestoreV2:tensors:33"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_33AssignVariableOpassignvariableop_33_t_kernelIdentity_33:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_34IdentityRestoreV2:tensors:34"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_34AssignVariableOpassignvariableop_34_t_biasIdentity_34:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_35IdentityRestoreV2:tensors:35"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_35AssignVariableOpassignvariableop_35_p_kernelIdentity_35:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_36IdentityRestoreV2:tensors:36"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_36AssignVariableOpassignvariableop_36_p_biasIdentity_36:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_37IdentityRestoreV2:tensors:37"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_37AssignVariableOpassignvariableop_37_phi_kernelIdentity_37:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_38IdentityRestoreV2:tensors:38"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_38AssignVariableOpassignvariableop_38_phi_biasIdentity_38:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_39IdentityRestoreV2:tensors:39"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_39AssignVariableOpassignvariableop_39_c_kernelIdentity_39:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_40IdentityRestoreV2:tensors:40"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_40AssignVariableOpassignvariableop_40_c_biasIdentity_40:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_41IdentityRestoreV2:tensors:41"/device:CPU:0*
T0	*
_output_shapes
:�
AssignVariableOp_41AssignVariableOpassignvariableop_41_iterationIdentity_41:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0	_
Identity_42IdentityRestoreV2:tensors:42"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_42AssignVariableOp!assignvariableop_42_learning_rateIdentity_42:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_43IdentityRestoreV2:tensors:43"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_43AssignVariableOp$assignvariableop_43_adam_m_l2_kernelIdentity_43:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_44IdentityRestoreV2:tensors:44"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_44AssignVariableOp$assignvariableop_44_adam_v_l2_kernelIdentity_44:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_45IdentityRestoreV2:tensors:45"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_45AssignVariableOp"assignvariableop_45_adam_m_l2_biasIdentity_45:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_46IdentityRestoreV2:tensors:46"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_46AssignVariableOp"assignvariableop_46_adam_v_l2_biasIdentity_46:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_47IdentityRestoreV2:tensors:47"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_47AssignVariableOp$assignvariableop_47_adam_m_l3_kernelIdentity_47:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_48IdentityRestoreV2:tensors:48"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_48AssignVariableOp$assignvariableop_48_adam_v_l3_kernelIdentity_48:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_49IdentityRestoreV2:tensors:49"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_49AssignVariableOp"assignvariableop_49_adam_m_l3_biasIdentity_49:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_50IdentityRestoreV2:tensors:50"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_50AssignVariableOp"assignvariableop_50_adam_v_l3_biasIdentity_50:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_51IdentityRestoreV2:tensors:51"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_51AssignVariableOp#assignvariableop_51_adam_m_x_kernelIdentity_51:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_52IdentityRestoreV2:tensors:52"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_52AssignVariableOp#assignvariableop_52_adam_v_x_kernelIdentity_52:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_53IdentityRestoreV2:tensors:53"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_53AssignVariableOp!assignvariableop_53_adam_m_x_biasIdentity_53:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_54IdentityRestoreV2:tensors:54"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_54AssignVariableOp!assignvariableop_54_adam_v_x_biasIdentity_54:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_55IdentityRestoreV2:tensors:55"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_55AssignVariableOp#assignvariableop_55_adam_m_1_kernelIdentity_55:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_56IdentityRestoreV2:tensors:56"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_56AssignVariableOp#assignvariableop_56_adam_v_1_kernelIdentity_56:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_57IdentityRestoreV2:tensors:57"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_57AssignVariableOp!assignvariableop_57_adam_m_1_biasIdentity_57:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_58IdentityRestoreV2:tensors:58"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_58AssignVariableOp!assignvariableop_58_adam_v_1_biasIdentity_58:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_59IdentityRestoreV2:tensors:59"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_59AssignVariableOp#assignvariableop_59_adam_m_4_kernelIdentity_59:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_60IdentityRestoreV2:tensors:60"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_60AssignVariableOp#assignvariableop_60_adam_v_4_kernelIdentity_60:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_61IdentityRestoreV2:tensors:61"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_61AssignVariableOp!assignvariableop_61_adam_m_4_biasIdentity_61:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_62IdentityRestoreV2:tensors:62"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_62AssignVariableOp!assignvariableop_62_adam_v_4_biasIdentity_62:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_63IdentityRestoreV2:tensors:63"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_63AssignVariableOp#assignvariableop_63_adam_m_7_kernelIdentity_63:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_64IdentityRestoreV2:tensors:64"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_64AssignVariableOp#assignvariableop_64_adam_v_7_kernelIdentity_64:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_65IdentityRestoreV2:tensors:65"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_65AssignVariableOp!assignvariableop_65_adam_m_7_biasIdentity_65:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_66IdentityRestoreV2:tensors:66"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_66AssignVariableOp!assignvariableop_66_adam_v_7_biasIdentity_66:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_67IdentityRestoreV2:tensors:67"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_67AssignVariableOp$assignvariableop_67_adam_m_10_kernelIdentity_67:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_68IdentityRestoreV2:tensors:68"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_68AssignVariableOp$assignvariableop_68_adam_v_10_kernelIdentity_68:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_69IdentityRestoreV2:tensors:69"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_69AssignVariableOp"assignvariableop_69_adam_m_10_biasIdentity_69:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_70IdentityRestoreV2:tensors:70"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_70AssignVariableOp"assignvariableop_70_adam_v_10_biasIdentity_70:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_71IdentityRestoreV2:tensors:71"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_71AssignVariableOp#assignvariableop_71_adam_m_2_kernelIdentity_71:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_72IdentityRestoreV2:tensors:72"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_72AssignVariableOp#assignvariableop_72_adam_v_2_kernelIdentity_72:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_73IdentityRestoreV2:tensors:73"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_73AssignVariableOp!assignvariableop_73_adam_m_2_biasIdentity_73:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_74IdentityRestoreV2:tensors:74"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_74AssignVariableOp!assignvariableop_74_adam_v_2_biasIdentity_74:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_75IdentityRestoreV2:tensors:75"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_75AssignVariableOp#assignvariableop_75_adam_m_5_kernelIdentity_75:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_76IdentityRestoreV2:tensors:76"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_76AssignVariableOp#assignvariableop_76_adam_v_5_kernelIdentity_76:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_77IdentityRestoreV2:tensors:77"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_77AssignVariableOp!assignvariableop_77_adam_m_5_biasIdentity_77:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_78IdentityRestoreV2:tensors:78"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_78AssignVariableOp!assignvariableop_78_adam_v_5_biasIdentity_78:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_79IdentityRestoreV2:tensors:79"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_79AssignVariableOp#assignvariableop_79_adam_m_8_kernelIdentity_79:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_80IdentityRestoreV2:tensors:80"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_80AssignVariableOp#assignvariableop_80_adam_v_8_kernelIdentity_80:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_81IdentityRestoreV2:tensors:81"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_81AssignVariableOp!assignvariableop_81_adam_m_8_biasIdentity_81:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_82IdentityRestoreV2:tensors:82"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_82AssignVariableOp!assignvariableop_82_adam_v_8_biasIdentity_82:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_83IdentityRestoreV2:tensors:83"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_83AssignVariableOp$assignvariableop_83_adam_m_11_kernelIdentity_83:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_84IdentityRestoreV2:tensors:84"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_84AssignVariableOp$assignvariableop_84_adam_v_11_kernelIdentity_84:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_85IdentityRestoreV2:tensors:85"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_85AssignVariableOp"assignvariableop_85_adam_m_11_biasIdentity_85:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_86IdentityRestoreV2:tensors:86"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_86AssignVariableOp"assignvariableop_86_adam_v_11_biasIdentity_86:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_87IdentityRestoreV2:tensors:87"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_87AssignVariableOp#assignvariableop_87_adam_m_3_kernelIdentity_87:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_88IdentityRestoreV2:tensors:88"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_88AssignVariableOp#assignvariableop_88_adam_v_3_kernelIdentity_88:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_89IdentityRestoreV2:tensors:89"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_89AssignVariableOp!assignvariableop_89_adam_m_3_biasIdentity_89:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_90IdentityRestoreV2:tensors:90"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_90AssignVariableOp!assignvariableop_90_adam_v_3_biasIdentity_90:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_91IdentityRestoreV2:tensors:91"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_91AssignVariableOp#assignvariableop_91_adam_m_6_kernelIdentity_91:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_92IdentityRestoreV2:tensors:92"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_92AssignVariableOp#assignvariableop_92_adam_v_6_kernelIdentity_92:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_93IdentityRestoreV2:tensors:93"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_93AssignVariableOp!assignvariableop_93_adam_m_6_biasIdentity_93:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_94IdentityRestoreV2:tensors:94"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_94AssignVariableOp!assignvariableop_94_adam_v_6_biasIdentity_94:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_95IdentityRestoreV2:tensors:95"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_95AssignVariableOp#assignvariableop_95_adam_m_9_kernelIdentity_95:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_96IdentityRestoreV2:tensors:96"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_96AssignVariableOp#assignvariableop_96_adam_v_9_kernelIdentity_96:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_97IdentityRestoreV2:tensors:97"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_97AssignVariableOp!assignvariableop_97_adam_m_9_biasIdentity_97:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_98IdentityRestoreV2:tensors:98"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_98AssignVariableOp!assignvariableop_98_adam_v_9_biasIdentity_98:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0_
Identity_99IdentityRestoreV2:tensors:99"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_99AssignVariableOp$assignvariableop_99_adam_m_12_kernelIdentity_99:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0a
Identity_100IdentityRestoreV2:tensors:100"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_100AssignVariableOp%assignvariableop_100_adam_v_12_kernelIdentity_100:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0a
Identity_101IdentityRestoreV2:tensors:101"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_101AssignVariableOp#assignvariableop_101_adam_m_12_biasIdentity_101:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0a
Identity_102IdentityRestoreV2:tensors:102"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_102AssignVariableOp#assignvariableop_102_adam_v_12_biasIdentity_102:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0a
Identity_103IdentityRestoreV2:tensors:103"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_103AssignVariableOp$assignvariableop_103_adam_m_t_kernelIdentity_103:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0a
Identity_104IdentityRestoreV2:tensors:104"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_104AssignVariableOp$assignvariableop_104_adam_v_t_kernelIdentity_104:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0a
Identity_105IdentityRestoreV2:tensors:105"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_105AssignVariableOp"assignvariableop_105_adam_m_t_biasIdentity_105:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0a
Identity_106IdentityRestoreV2:tensors:106"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_106AssignVariableOp"assignvariableop_106_adam_v_t_biasIdentity_106:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0a
Identity_107IdentityRestoreV2:tensors:107"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_107AssignVariableOp$assignvariableop_107_adam_m_p_kernelIdentity_107:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0a
Identity_108IdentityRestoreV2:tensors:108"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_108AssignVariableOp$assignvariableop_108_adam_v_p_kernelIdentity_108:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0a
Identity_109IdentityRestoreV2:tensors:109"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_109AssignVariableOp"assignvariableop_109_adam_m_p_biasIdentity_109:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0a
Identity_110IdentityRestoreV2:tensors:110"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_110AssignVariableOp"assignvariableop_110_adam_v_p_biasIdentity_110:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0a
Identity_111IdentityRestoreV2:tensors:111"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_111AssignVariableOp&assignvariableop_111_adam_m_phi_kernelIdentity_111:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0a
Identity_112IdentityRestoreV2:tensors:112"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_112AssignVariableOp&assignvariableop_112_adam_v_phi_kernelIdentity_112:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0a
Identity_113IdentityRestoreV2:tensors:113"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_113AssignVariableOp$assignvariableop_113_adam_m_phi_biasIdentity_113:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0a
Identity_114IdentityRestoreV2:tensors:114"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_114AssignVariableOp$assignvariableop_114_adam_v_phi_biasIdentity_114:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0a
Identity_115IdentityRestoreV2:tensors:115"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_115AssignVariableOp$assignvariableop_115_adam_m_c_kernelIdentity_115:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0a
Identity_116IdentityRestoreV2:tensors:116"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_116AssignVariableOp$assignvariableop_116_adam_v_c_kernelIdentity_116:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0a
Identity_117IdentityRestoreV2:tensors:117"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_117AssignVariableOp"assignvariableop_117_adam_m_c_biasIdentity_117:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0a
Identity_118IdentityRestoreV2:tensors:118"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_118AssignVariableOp"assignvariableop_118_adam_v_c_biasIdentity_118:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0a
Identity_119IdentityRestoreV2:tensors:119"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_119AssignVariableOpassignvariableop_119_total_4Identity_119:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0a
Identity_120IdentityRestoreV2:tensors:120"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_120AssignVariableOpassignvariableop_120_count_4Identity_120:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0a
Identity_121IdentityRestoreV2:tensors:121"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_121AssignVariableOpassignvariableop_121_total_3Identity_121:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0a
Identity_122IdentityRestoreV2:tensors:122"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_122AssignVariableOpassignvariableop_122_count_3Identity_122:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0a
Identity_123IdentityRestoreV2:tensors:123"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_123AssignVariableOpassignvariableop_123_total_2Identity_123:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0a
Identity_124IdentityRestoreV2:tensors:124"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_124AssignVariableOpassignvariableop_124_count_2Identity_124:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0a
Identity_125IdentityRestoreV2:tensors:125"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_125AssignVariableOpassignvariableop_125_total_1Identity_125:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0a
Identity_126IdentityRestoreV2:tensors:126"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_126AssignVariableOpassignvariableop_126_count_1Identity_126:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0a
Identity_127IdentityRestoreV2:tensors:127"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_127AssignVariableOpassignvariableop_127_totalIdentity_127:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0a
Identity_128IdentityRestoreV2:tensors:128"/device:CPU:0*
T0*
_output_shapes
:�
AssignVariableOp_128AssignVariableOpassignvariableop_128_countIdentity_128:output:0"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *
dtype0Y
NoOpNoOp"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 �
Identity_129Identityfile_prefix^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_100^AssignVariableOp_101^AssignVariableOp_102^AssignVariableOp_103^AssignVariableOp_104^AssignVariableOp_105^AssignVariableOp_106^AssignVariableOp_107^AssignVariableOp_108^AssignVariableOp_109^AssignVariableOp_11^AssignVariableOp_110^AssignVariableOp_111^AssignVariableOp_112^AssignVariableOp_113^AssignVariableOp_114^AssignVariableOp_115^AssignVariableOp_116^AssignVariableOp_117^AssignVariableOp_118^AssignVariableOp_119^AssignVariableOp_12^AssignVariableOp_120^AssignVariableOp_121^AssignVariableOp_122^AssignVariableOp_123^AssignVariableOp_124^AssignVariableOp_125^AssignVariableOp_126^AssignVariableOp_127^AssignVariableOp_128^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_25^AssignVariableOp_26^AssignVariableOp_27^AssignVariableOp_28^AssignVariableOp_29^AssignVariableOp_3^AssignVariableOp_30^AssignVariableOp_31^AssignVariableOp_32^AssignVariableOp_33^AssignVariableOp_34^AssignVariableOp_35^AssignVariableOp_36^AssignVariableOp_37^AssignVariableOp_38^AssignVariableOp_39^AssignVariableOp_4^AssignVariableOp_40^AssignVariableOp_41^AssignVariableOp_42^AssignVariableOp_43^AssignVariableOp_44^AssignVariableOp_45^AssignVariableOp_46^AssignVariableOp_47^AssignVariableOp_48^AssignVariableOp_49^AssignVariableOp_5^AssignVariableOp_50^AssignVariableOp_51^AssignVariableOp_52^AssignVariableOp_53^AssignVariableOp_54^AssignVariableOp_55^AssignVariableOp_56^AssignVariableOp_57^AssignVariableOp_58^AssignVariableOp_59^AssignVariableOp_6^AssignVariableOp_60^AssignVariableOp_61^AssignVariableOp_62^AssignVariableOp_63^AssignVariableOp_64^AssignVariableOp_65^AssignVariableOp_66^AssignVariableOp_67^AssignVariableOp_68^AssignVariableOp_69^AssignVariableOp_7^AssignVariableOp_70^AssignVariableOp_71^AssignVariableOp_72^AssignVariableOp_73^AssignVariableOp_74^AssignVariableOp_75^AssignVariableOp_76^AssignVariableOp_77^AssignVariableOp_78^AssignVariableOp_79^AssignVariableOp_8^AssignVariableOp_80^AssignVariableOp_81^AssignVariableOp_82^AssignVariableOp_83^AssignVariableOp_84^AssignVariableOp_85^AssignVariableOp_86^AssignVariableOp_87^AssignVariableOp_88^AssignVariableOp_89^AssignVariableOp_9^AssignVariableOp_90^AssignVariableOp_91^AssignVariableOp_92^AssignVariableOp_93^AssignVariableOp_94^AssignVariableOp_95^AssignVariableOp_96^AssignVariableOp_97^AssignVariableOp_98^AssignVariableOp_99^NoOp"/device:CPU:0*
T0*
_output_shapes
: Y
Identity_130IdentityIdentity_129:output:0^NoOp_1*
T0*
_output_shapes
: �
NoOp_1NoOp^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_100^AssignVariableOp_101^AssignVariableOp_102^AssignVariableOp_103^AssignVariableOp_104^AssignVariableOp_105^AssignVariableOp_106^AssignVariableOp_107^AssignVariableOp_108^AssignVariableOp_109^AssignVariableOp_11^AssignVariableOp_110^AssignVariableOp_111^AssignVariableOp_112^AssignVariableOp_113^AssignVariableOp_114^AssignVariableOp_115^AssignVariableOp_116^AssignVariableOp_117^AssignVariableOp_118^AssignVariableOp_119^AssignVariableOp_12^AssignVariableOp_120^AssignVariableOp_121^AssignVariableOp_122^AssignVariableOp_123^AssignVariableOp_124^AssignVariableOp_125^AssignVariableOp_126^AssignVariableOp_127^AssignVariableOp_128^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_20^AssignVariableOp_21^AssignVariableOp_22^AssignVariableOp_23^AssignVariableOp_24^AssignVariableOp_25^AssignVariableOp_26^AssignVariableOp_27^AssignVariableOp_28^AssignVariableOp_29^AssignVariableOp_3^AssignVariableOp_30^AssignVariableOp_31^AssignVariableOp_32^AssignVariableOp_33^AssignVariableOp_34^AssignVariableOp_35^AssignVariableOp_36^AssignVariableOp_37^AssignVariableOp_38^AssignVariableOp_39^AssignVariableOp_4^AssignVariableOp_40^AssignVariableOp_41^AssignVariableOp_42^AssignVariableOp_43^AssignVariableOp_44^AssignVariableOp_45^AssignVariableOp_46^AssignVariableOp_47^AssignVariableOp_48^AssignVariableOp_49^AssignVariableOp_5^AssignVariableOp_50^AssignVariableOp_51^AssignVariableOp_52^AssignVariableOp_53^AssignVariableOp_54^AssignVariableOp_55^AssignVariableOp_56^AssignVariableOp_57^AssignVariableOp_58^AssignVariableOp_59^AssignVariableOp_6^AssignVariableOp_60^AssignVariableOp_61^AssignVariableOp_62^AssignVariableOp_63^AssignVariableOp_64^AssignVariableOp_65^AssignVariableOp_66^AssignVariableOp_67^AssignVariableOp_68^AssignVariableOp_69^AssignVariableOp_7^AssignVariableOp_70^AssignVariableOp_71^AssignVariableOp_72^AssignVariableOp_73^AssignVariableOp_74^AssignVariableOp_75^AssignVariableOp_76^AssignVariableOp_77^AssignVariableOp_78^AssignVariableOp_79^AssignVariableOp_8^AssignVariableOp_80^AssignVariableOp_81^AssignVariableOp_82^AssignVariableOp_83^AssignVariableOp_84^AssignVariableOp_85^AssignVariableOp_86^AssignVariableOp_87^AssignVariableOp_88^AssignVariableOp_89^AssignVariableOp_9^AssignVariableOp_90^AssignVariableOp_91^AssignVariableOp_92^AssignVariableOp_93^AssignVariableOp_94^AssignVariableOp_95^AssignVariableOp_96^AssignVariableOp_97^AssignVariableOp_98^AssignVariableOp_99*"
_acd_function_control_output(*
_output_shapes
 "%
identity_130Identity_130:output:0*�
_input_shapes�
�: : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : 2$
AssignVariableOpAssignVariableOp2(
AssignVariableOp_1AssignVariableOp_12*
AssignVariableOp_10AssignVariableOp_102,
AssignVariableOp_100AssignVariableOp_1002,
AssignVariableOp_101AssignVariableOp_1012,
AssignVariableOp_102AssignVariableOp_1022,
AssignVariableOp_103AssignVariableOp_1032,
AssignVariableOp_104AssignVariableOp_1042,
AssignVariableOp_105AssignVariableOp_1052,
AssignVariableOp_106AssignVariableOp_1062,
AssignVariableOp_107AssignVariableOp_1072,
AssignVariableOp_108AssignVariableOp_1082,
AssignVariableOp_109AssignVariableOp_1092*
AssignVariableOp_11AssignVariableOp_112,
AssignVariableOp_110AssignVariableOp_1102,
AssignVariableOp_111AssignVariableOp_1112,
AssignVariableOp_112AssignVariableOp_1122,
AssignVariableOp_113AssignVariableOp_1132,
AssignVariableOp_114AssignVariableOp_1142,
AssignVariableOp_115AssignVariableOp_1152,
AssignVariableOp_116AssignVariableOp_1162,
AssignVariableOp_117AssignVariableOp_1172,
AssignVariableOp_118AssignVariableOp_1182,
AssignVariableOp_119AssignVariableOp_1192*
AssignVariableOp_12AssignVariableOp_122,
AssignVariableOp_120AssignVariableOp_1202,
AssignVariableOp_121AssignVariableOp_1212,
AssignVariableOp_122AssignVariableOp_1222,
AssignVariableOp_123AssignVariableOp_1232,
AssignVariableOp_124AssignVariableOp_1242,
AssignVariableOp_125AssignVariableOp_1252,
AssignVariableOp_126AssignVariableOp_1262,
AssignVariableOp_127AssignVariableOp_1272,
AssignVariableOp_128AssignVariableOp_1282*
AssignVariableOp_13AssignVariableOp_132*
AssignVariableOp_14AssignVariableOp_142*
AssignVariableOp_15AssignVariableOp_152*
AssignVariableOp_16AssignVariableOp_162*
AssignVariableOp_17AssignVariableOp_172*
AssignVariableOp_18AssignVariableOp_182*
AssignVariableOp_19AssignVariableOp_192(
AssignVariableOp_2AssignVariableOp_22*
AssignVariableOp_20AssignVariableOp_202*
AssignVariableOp_21AssignVariableOp_212*
AssignVariableOp_22AssignVariableOp_222*
AssignVariableOp_23AssignVariableOp_232*
AssignVariableOp_24AssignVariableOp_242*
AssignVariableOp_25AssignVariableOp_252*
AssignVariableOp_26AssignVariableOp_262*
AssignVariableOp_27AssignVariableOp_272*
AssignVariableOp_28AssignVariableOp_282*
AssignVariableOp_29AssignVariableOp_292(
AssignVariableOp_3AssignVariableOp_32*
AssignVariableOp_30AssignVariableOp_302*
AssignVariableOp_31AssignVariableOp_312*
AssignVariableOp_32AssignVariableOp_322*
AssignVariableOp_33AssignVariableOp_332*
AssignVariableOp_34AssignVariableOp_342*
AssignVariableOp_35AssignVariableOp_352*
AssignVariableOp_36AssignVariableOp_362*
AssignVariableOp_37AssignVariableOp_372*
AssignVariableOp_38AssignVariableOp_382*
AssignVariableOp_39AssignVariableOp_392(
AssignVariableOp_4AssignVariableOp_42*
AssignVariableOp_40AssignVariableOp_402*
AssignVariableOp_41AssignVariableOp_412*
AssignVariableOp_42AssignVariableOp_422*
AssignVariableOp_43AssignVariableOp_432*
AssignVariableOp_44AssignVariableOp_442*
AssignVariableOp_45AssignVariableOp_452*
AssignVariableOp_46AssignVariableOp_462*
AssignVariableOp_47AssignVariableOp_472*
AssignVariableOp_48AssignVariableOp_482*
AssignVariableOp_49AssignVariableOp_492(
AssignVariableOp_5AssignVariableOp_52*
AssignVariableOp_50AssignVariableOp_502*
AssignVariableOp_51AssignVariableOp_512*
AssignVariableOp_52AssignVariableOp_522*
AssignVariableOp_53AssignVariableOp_532*
AssignVariableOp_54AssignVariableOp_542*
AssignVariableOp_55AssignVariableOp_552*
AssignVariableOp_56AssignVariableOp_562*
AssignVariableOp_57AssignVariableOp_572*
AssignVariableOp_58AssignVariableOp_582*
AssignVariableOp_59AssignVariableOp_592(
AssignVariableOp_6AssignVariableOp_62*
AssignVariableOp_60AssignVariableOp_602*
AssignVariableOp_61AssignVariableOp_612*
AssignVariableOp_62AssignVariableOp_622*
AssignVariableOp_63AssignVariableOp_632*
AssignVariableOp_64AssignVariableOp_642*
AssignVariableOp_65AssignVariableOp_652*
AssignVariableOp_66AssignVariableOp_662*
AssignVariableOp_67AssignVariableOp_672*
AssignVariableOp_68AssignVariableOp_682*
AssignVariableOp_69AssignVariableOp_692(
AssignVariableOp_7AssignVariableOp_72*
AssignVariableOp_70AssignVariableOp_702*
AssignVariableOp_71AssignVariableOp_712*
AssignVariableOp_72AssignVariableOp_722*
AssignVariableOp_73AssignVariableOp_732*
AssignVariableOp_74AssignVariableOp_742*
AssignVariableOp_75AssignVariableOp_752*
AssignVariableOp_76AssignVariableOp_762*
AssignVariableOp_77AssignVariableOp_772*
AssignVariableOp_78AssignVariableOp_782*
AssignVariableOp_79AssignVariableOp_792(
AssignVariableOp_8AssignVariableOp_82*
AssignVariableOp_80AssignVariableOp_802*
AssignVariableOp_81AssignVariableOp_812*
AssignVariableOp_82AssignVariableOp_822*
AssignVariableOp_83AssignVariableOp_832*
AssignVariableOp_84AssignVariableOp_842*
AssignVariableOp_85AssignVariableOp_852*
AssignVariableOp_86AssignVariableOp_862*
AssignVariableOp_87AssignVariableOp_872*
AssignVariableOp_88AssignVariableOp_882*
AssignVariableOp_89AssignVariableOp_892(
AssignVariableOp_9AssignVariableOp_92*
AssignVariableOp_90AssignVariableOp_902*
AssignVariableOp_91AssignVariableOp_912*
AssignVariableOp_92AssignVariableOp_922*
AssignVariableOp_93AssignVariableOp_932*
AssignVariableOp_94AssignVariableOp_942*
AssignVariableOp_95AssignVariableOp_952*
AssignVariableOp_96AssignVariableOp_962*
AssignVariableOp_97AssignVariableOp_972*
AssignVariableOp_98AssignVariableOp_982*
AssignVariableOp_99AssignVariableOp_99:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix
�

�
>__inference_2_layer_call_and_return_conditional_losses_2218860

inputs0
matmul_readvariableop_resource:@@-
biasadd_readvariableop_resource:@
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@@*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������@a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:���������@w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�

�
>__inference_8_layer_call_and_return_conditional_losses_2218900

inputs0
matmul_readvariableop_resource:@@-
biasadd_readvariableop_resource:@
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@@*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������@a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:���������@w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�

�
>__inference_3_layer_call_and_return_conditional_losses_2217127

inputs0
matmul_readvariableop_resource:@ -
biasadd_readvariableop_resource: 
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@ *
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:��������� a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:��������� w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�
�
$__inference_l3_layer_call_fn_2218729

inputs
unknown:@@
	unknown_0:@
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *H
fCRA
?__inference_l3_layer_call_and_return_conditional_losses_2216906o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������@`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�

�
>__inference_3_layer_call_and_return_conditional_losses_2218940

inputs0
matmul_readvariableop_resource:@ -
biasadd_readvariableop_resource: 
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@ *
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:��������� a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:��������� w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�
�
#__inference_8_layer_call_fn_2218889

inputs
unknown:@@
	unknown_0:@
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_8_layer_call_and_return_conditional_losses_2217025o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������@`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�
�
#__inference_6_layer_call_fn_2218949

inputs
unknown:@ 
	unknown_0: 
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:��������� *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_6_layer_call_and_return_conditional_losses_2217110o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:��������� `
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�

�
>__inference_P_layer_call_and_return_conditional_losses_2219040

inputs0
matmul_readvariableop_resource: -
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

: *
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:���������w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:��������� : : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:��������� 
 
_user_specified_nameinputs
�

�
@__inference_phi_layer_call_and_return_conditional_losses_2219060

inputs0
matmul_readvariableop_resource: -
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

: *
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:���������w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:��������� : : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:��������� 
 
_user_specified_nameinputs
�'
�
__inference_adapt_step_2218224
iterator%
add_readvariableop_resource:	 %
readvariableop_resource:'
readvariableop_2_resource:��AssignVariableOp�AssignVariableOp_1�AssignVariableOp_2�IteratorGetNext�ReadVariableOp�ReadVariableOp_1�ReadVariableOp_2�add/ReadVariableOp�
IteratorGetNextIteratorGetNextiterator*
_class
loc:@iterator*'
_output_shapes
:���������*&
output_shapes
:���������*
output_types
2h
moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: �
moments/meanMeanIteratorGetNext:components:0'moments/mean/reduction_indices:output:0*
T0*
_output_shapes

:*
	keep_dims(d
moments/StopGradientStopGradientmoments/mean:output:0*
T0*
_output_shapes

:�
moments/SquaredDifferenceSquaredDifferenceIteratorGetNext:components:0moments/StopGradient:output:0*
T0*'
_output_shapes
:���������l
"moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: �
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
T0*
_output_shapes

:*
	keep_dims(m
moments/SqueezeSqueezemoments/mean:output:0*
T0*
_output_shapes
:*
squeeze_dims
 s
moments/Squeeze_1Squeezemoments/variance:output:0*
T0*
_output_shapes
:*
squeeze_dims
 a
ShapeShapeIteratorGetNext:components:0*
T0*
_output_shapes
:*
out_type0	Z
GatherV2/indicesConst*
_output_shapes
:*
dtype0*
valueB: O
GatherV2/axisConst*
_output_shapes
: *
dtype0*
value	B : �
GatherV2GatherV2Shape:output:0GatherV2/indices:output:0GatherV2/axis:output:0*
Taxis0*
Tindices0*
Tparams0	*
_output_shapes
:O
ConstConst*
_output_shapes
:*
dtype0*
valueB: P
ProdProdGatherV2:output:0Const:output:0*
T0	*
_output_shapes
: f
add/ReadVariableOpReadVariableOpadd_readvariableop_resource*
_output_shapes
: *
dtype0	X
addAddV2Prod:output:0add/ReadVariableOp:value:0*
T0	*
_output_shapes
: K
CastCastProd:output:0*

DstT0*

SrcT0	*
_output_shapes
: G
Cast_1Castadd:z:0*

DstT0*

SrcT0	*
_output_shapes
: I
truedivRealDivCast:y:0
Cast_1:y:0*
T0*
_output_shapes
: J
sub/xConst*
_output_shapes
: *
dtype0*
valueB
 *  �?H
subSubsub/x:output:0truediv:z:0*
T0*
_output_shapes
: b
ReadVariableOpReadVariableOpreadvariableop_resource*
_output_shapes
:*
dtype0P
mulMulReadVariableOp:value:0sub:z:0*
T0*
_output_shapes
:X
mul_1Mulmoments/Squeeze:output:0truediv:z:0*
T0*
_output_shapes
:G
add_1AddV2mul:z:0	mul_1:z:0*
T0*
_output_shapes
:d
ReadVariableOp_1ReadVariableOpreadvariableop_resource*
_output_shapes
:*
dtype0V
sub_1SubReadVariableOp_1:value:0	add_1:z:0*
T0*
_output_shapes
:J
pow/yConst*
_output_shapes
: *
dtype0*
valueB
 *   @J
powPow	sub_1:z:0pow/y:output:0*
T0*
_output_shapes
:f
ReadVariableOp_2ReadVariableOpreadvariableop_2_resource*
_output_shapes
:*
dtype0V
add_2AddV2ReadVariableOp_2:value:0pow:z:0*
T0*
_output_shapes
:E
mul_2Mul	add_2:z:0sub:z:0*
T0*
_output_shapes
:V
sub_2Submoments/Squeeze:output:0	add_1:z:0*
T0*
_output_shapes
:L
pow_1/yConst*
_output_shapes
: *
dtype0*
valueB
 *   @N
pow_1Pow	sub_2:z:0pow_1/y:output:0*
T0*
_output_shapes
:Z
add_3AddV2moments/Squeeze_1:output:0	pow_1:z:0*
T0*
_output_shapes
:I
mul_3Mul	add_3:z:0truediv:z:0*
T0*
_output_shapes
:I
add_4AddV2	mul_2:z:0	mul_3:z:0*
T0*
_output_shapes
:�
AssignVariableOpAssignVariableOpreadvariableop_resource	add_1:z:0^ReadVariableOp^ReadVariableOp_1*
_output_shapes
 *
dtype0*
validate_shape(�
AssignVariableOp_1AssignVariableOpreadvariableop_2_resource	add_4:z:0^ReadVariableOp_2*
_output_shapes
 *
dtype0*
validate_shape(�
AssignVariableOp_2AssignVariableOpadd_readvariableop_resourceadd:z:0^add/ReadVariableOp*
_output_shapes
 *
dtype0	*
validate_shape(*(
_construction_contextkEagerRuntime*
_input_shapes

: : : : 2$
AssignVariableOpAssignVariableOp2(
AssignVariableOp_1AssignVariableOp_12(
AssignVariableOp_2AssignVariableOp_22"
IteratorGetNextIteratorGetNext2 
ReadVariableOpReadVariableOp2$
ReadVariableOp_1ReadVariableOp_12$
ReadVariableOp_2ReadVariableOp_22(
add/ReadVariableOpadd/ReadVariableOp:( $
"
_user_specified_name
iterator
�^
�
B__inference_model_layer_call_and_return_conditional_losses_2217975	
input
normalization_sub_y
normalization_sqrt_x

l2_2217876:@

l2_2217878:@

l3_2217881:@@

l3_2217883:@
	x_2217886:@@
	x_2217888:@
unknown:@@
	unknown_0:@
	unknown_1:@@
	unknown_2:@
	unknown_3:@@
	unknown_4:@
	unknown_5:@@
	unknown_6:@
	unknown_7:@@
	unknown_8:@
	unknown_9:@@

unknown_10:@

unknown_11:@@

unknown_12:@

unknown_13:@@

unknown_14:@

unknown_15:@ 

unknown_16: 

unknown_17:@ 

unknown_18: 

unknown_19:@ 

unknown_20: 

unknown_21:@ 

unknown_22: 
	c_2217951: 
	c_2217953:
phi_2217956: 
phi_2217958:
	p_2217961: 
	p_2217963:
	t_2217966: 
	t_2217968:
identity

identity_1

identity_2

identity_3��1/StatefulPartitionedCall�10/StatefulPartitionedCall�11/StatefulPartitionedCall�12/StatefulPartitionedCall�2/StatefulPartitionedCall�3/StatefulPartitionedCall�4/StatefulPartitionedCall�5/StatefulPartitionedCall�6/StatefulPartitionedCall�7/StatefulPartitionedCall�8/StatefulPartitionedCall�9/StatefulPartitionedCall�P/StatefulPartitionedCall�T/StatefulPartitionedCall�c/StatefulPartitionedCall�l2/StatefulPartitionedCall�l3/StatefulPartitionedCall�phi/StatefulPartitionedCall�x/StatefulPartitionedCallf
normalization/subSubinputnormalization_sub_y*
T0*'
_output_shapes
:���������Y
normalization/SqrtSqrtnormalization_sqrt_x*
T0*
_output_shapes

:\
normalization/Maximum/yConst*
_output_shapes
: *
dtype0*
valueB
 *���3�
normalization/MaximumMaximumnormalization/Sqrt:y:0 normalization/Maximum/y:output:0*
T0*
_output_shapes

:�
normalization/truedivRealDivnormalization/sub:z:0normalization/Maximum:z:0*
T0*'
_output_shapes
:����������
l2/StatefulPartitionedCallStatefulPartitionedCallnormalization/truediv:z:0
l2_2217876
l2_2217878*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *H
fCRA
?__inference_l2_layer_call_and_return_conditional_losses_2216889�
l3/StatefulPartitionedCallStatefulPartitionedCall#l2/StatefulPartitionedCall:output:0
l3_2217881
l3_2217883*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *H
fCRA
?__inference_l3_layer_call_and_return_conditional_losses_2216906�
x/StatefulPartitionedCallStatefulPartitionedCall#l3/StatefulPartitionedCall:output:0	x_2217886	x_2217888*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_x_layer_call_and_return_conditional_losses_2216923�
10/StatefulPartitionedCallStatefulPartitionedCall"x/StatefulPartitionedCall:output:0unknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *H
fCRA
?__inference_10_layer_call_and_return_conditional_losses_2216940�
7/StatefulPartitionedCallStatefulPartitionedCall"x/StatefulPartitionedCall:output:0	unknown_1	unknown_2*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_7_layer_call_and_return_conditional_losses_2216957�
4/StatefulPartitionedCallStatefulPartitionedCall"x/StatefulPartitionedCall:output:0	unknown_3	unknown_4*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_4_layer_call_and_return_conditional_losses_2216974�
1/StatefulPartitionedCallStatefulPartitionedCall"x/StatefulPartitionedCall:output:0	unknown_5	unknown_6*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_1_layer_call_and_return_conditional_losses_2216991�
11/StatefulPartitionedCallStatefulPartitionedCall#10/StatefulPartitionedCall:output:0	unknown_7	unknown_8*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *H
fCRA
?__inference_11_layer_call_and_return_conditional_losses_2217008�
8/StatefulPartitionedCallStatefulPartitionedCall"7/StatefulPartitionedCall:output:0	unknown_9
unknown_10*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_8_layer_call_and_return_conditional_losses_2217025�
5/StatefulPartitionedCallStatefulPartitionedCall"4/StatefulPartitionedCall:output:0
unknown_11
unknown_12*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_5_layer_call_and_return_conditional_losses_2217042�
2/StatefulPartitionedCallStatefulPartitionedCall"1/StatefulPartitionedCall:output:0
unknown_13
unknown_14*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_2_layer_call_and_return_conditional_losses_2217059�
12/StatefulPartitionedCallStatefulPartitionedCall#11/StatefulPartitionedCall:output:0
unknown_15
unknown_16*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:��������� *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *H
fCRA
?__inference_12_layer_call_and_return_conditional_losses_2217076�
9/StatefulPartitionedCallStatefulPartitionedCall"8/StatefulPartitionedCall:output:0
unknown_17
unknown_18*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:��������� *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_9_layer_call_and_return_conditional_losses_2217093�
6/StatefulPartitionedCallStatefulPartitionedCall"5/StatefulPartitionedCall:output:0
unknown_19
unknown_20*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:��������� *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_6_layer_call_and_return_conditional_losses_2217110�
3/StatefulPartitionedCallStatefulPartitionedCall"2/StatefulPartitionedCall:output:0
unknown_21
unknown_22*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:��������� *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_3_layer_call_and_return_conditional_losses_2217127�
c/StatefulPartitionedCallStatefulPartitionedCall#12/StatefulPartitionedCall:output:0	c_2217951	c_2217953*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_c_layer_call_and_return_conditional_losses_2217144�
phi/StatefulPartitionedCallStatefulPartitionedCall"9/StatefulPartitionedCall:output:0phi_2217956phi_2217958*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *I
fDRB
@__inference_phi_layer_call_and_return_conditional_losses_2217161�
P/StatefulPartitionedCallStatefulPartitionedCall"6/StatefulPartitionedCall:output:0	p_2217961	p_2217963*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_P_layer_call_and_return_conditional_losses_2217178�
T/StatefulPartitionedCallStatefulPartitionedCall"3/StatefulPartitionedCall:output:0	t_2217966	t_2217968*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_T_layer_call_and_return_conditional_losses_2217195q
IdentityIdentity"T/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������s

Identity_1Identity"P/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������u

Identity_2Identity$phi/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������s

Identity_3Identity"c/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp^1/StatefulPartitionedCall^10/StatefulPartitionedCall^11/StatefulPartitionedCall^12/StatefulPartitionedCall^2/StatefulPartitionedCall^3/StatefulPartitionedCall^4/StatefulPartitionedCall^5/StatefulPartitionedCall^6/StatefulPartitionedCall^7/StatefulPartitionedCall^8/StatefulPartitionedCall^9/StatefulPartitionedCall^P/StatefulPartitionedCall^T/StatefulPartitionedCall^c/StatefulPartitionedCall^l2/StatefulPartitionedCall^l3/StatefulPartitionedCall^phi/StatefulPartitionedCall^x/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0"!

identity_3Identity_3:output:0*(
_construction_contextkEagerRuntime*�
_input_shapesu
s:���������::: : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : 26
1/StatefulPartitionedCall1/StatefulPartitionedCall28
10/StatefulPartitionedCall10/StatefulPartitionedCall28
11/StatefulPartitionedCall11/StatefulPartitionedCall28
12/StatefulPartitionedCall12/StatefulPartitionedCall26
2/StatefulPartitionedCall2/StatefulPartitionedCall26
3/StatefulPartitionedCall3/StatefulPartitionedCall26
4/StatefulPartitionedCall4/StatefulPartitionedCall26
5/StatefulPartitionedCall5/StatefulPartitionedCall26
6/StatefulPartitionedCall6/StatefulPartitionedCall26
7/StatefulPartitionedCall7/StatefulPartitionedCall26
8/StatefulPartitionedCall8/StatefulPartitionedCall26
9/StatefulPartitionedCall9/StatefulPartitionedCall26
P/StatefulPartitionedCallP/StatefulPartitionedCall26
T/StatefulPartitionedCallT/StatefulPartitionedCall26
c/StatefulPartitionedCallc/StatefulPartitionedCall28
l2/StatefulPartitionedCalll2/StatefulPartitionedCall28
l3/StatefulPartitionedCalll3/StatefulPartitionedCall2:
phi/StatefulPartitionedCallphi/StatefulPartitionedCall26
x/StatefulPartitionedCallx/StatefulPartitionedCall:N J
'
_output_shapes
:���������

_user_specified_nameinput:$ 

_output_shapes

::$ 

_output_shapes

:
�

�
?__inference_11_layer_call_and_return_conditional_losses_2217008

inputs0
matmul_readvariableop_resource:@@-
biasadd_readvariableop_resource:@
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@@*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������@a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:���������@w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
��
�-
 __inference__traced_save_2219495
file_prefix#
savev2_mean_read_readvariableop'
#savev2_variance_read_readvariableop&
"savev2_count_5_read_readvariableop	(
$savev2_l2_kernel_read_readvariableop&
"savev2_l2_bias_read_readvariableop(
$savev2_l3_kernel_read_readvariableop&
"savev2_l3_bias_read_readvariableop'
#savev2_x_kernel_read_readvariableop%
!savev2_x_bias_read_readvariableop'
#savev2_1_kernel_read_readvariableop%
!savev2_1_bias_read_readvariableop'
#savev2_4_kernel_read_readvariableop%
!savev2_4_bias_read_readvariableop'
#savev2_7_kernel_read_readvariableop%
!savev2_7_bias_read_readvariableop(
$savev2_10_kernel_read_readvariableop&
"savev2_10_bias_read_readvariableop'
#savev2_2_kernel_read_readvariableop%
!savev2_2_bias_read_readvariableop'
#savev2_5_kernel_read_readvariableop%
!savev2_5_bias_read_readvariableop'
#savev2_8_kernel_read_readvariableop%
!savev2_8_bias_read_readvariableop(
$savev2_11_kernel_read_readvariableop&
"savev2_11_bias_read_readvariableop'
#savev2_3_kernel_read_readvariableop%
!savev2_3_bias_read_readvariableop'
#savev2_6_kernel_read_readvariableop%
!savev2_6_bias_read_readvariableop'
#savev2_9_kernel_read_readvariableop%
!savev2_9_bias_read_readvariableop(
$savev2_12_kernel_read_readvariableop&
"savev2_12_bias_read_readvariableop'
#savev2_t_kernel_read_readvariableop%
!savev2_t_bias_read_readvariableop'
#savev2_p_kernel_read_readvariableop%
!savev2_p_bias_read_readvariableop)
%savev2_phi_kernel_read_readvariableop'
#savev2_phi_bias_read_readvariableop'
#savev2_c_kernel_read_readvariableop%
!savev2_c_bias_read_readvariableop(
$savev2_iteration_read_readvariableop	,
(savev2_learning_rate_read_readvariableop/
+savev2_adam_m_l2_kernel_read_readvariableop/
+savev2_adam_v_l2_kernel_read_readvariableop-
)savev2_adam_m_l2_bias_read_readvariableop-
)savev2_adam_v_l2_bias_read_readvariableop/
+savev2_adam_m_l3_kernel_read_readvariableop/
+savev2_adam_v_l3_kernel_read_readvariableop-
)savev2_adam_m_l3_bias_read_readvariableop-
)savev2_adam_v_l3_bias_read_readvariableop.
*savev2_adam_m_x_kernel_read_readvariableop.
*savev2_adam_v_x_kernel_read_readvariableop,
(savev2_adam_m_x_bias_read_readvariableop,
(savev2_adam_v_x_bias_read_readvariableop.
*savev2_adam_m_1_kernel_read_readvariableop.
*savev2_adam_v_1_kernel_read_readvariableop,
(savev2_adam_m_1_bias_read_readvariableop,
(savev2_adam_v_1_bias_read_readvariableop.
*savev2_adam_m_4_kernel_read_readvariableop.
*savev2_adam_v_4_kernel_read_readvariableop,
(savev2_adam_m_4_bias_read_readvariableop,
(savev2_adam_v_4_bias_read_readvariableop.
*savev2_adam_m_7_kernel_read_readvariableop.
*savev2_adam_v_7_kernel_read_readvariableop,
(savev2_adam_m_7_bias_read_readvariableop,
(savev2_adam_v_7_bias_read_readvariableop/
+savev2_adam_m_10_kernel_read_readvariableop/
+savev2_adam_v_10_kernel_read_readvariableop-
)savev2_adam_m_10_bias_read_readvariableop-
)savev2_adam_v_10_bias_read_readvariableop.
*savev2_adam_m_2_kernel_read_readvariableop.
*savev2_adam_v_2_kernel_read_readvariableop,
(savev2_adam_m_2_bias_read_readvariableop,
(savev2_adam_v_2_bias_read_readvariableop.
*savev2_adam_m_5_kernel_read_readvariableop.
*savev2_adam_v_5_kernel_read_readvariableop,
(savev2_adam_m_5_bias_read_readvariableop,
(savev2_adam_v_5_bias_read_readvariableop.
*savev2_adam_m_8_kernel_read_readvariableop.
*savev2_adam_v_8_kernel_read_readvariableop,
(savev2_adam_m_8_bias_read_readvariableop,
(savev2_adam_v_8_bias_read_readvariableop/
+savev2_adam_m_11_kernel_read_readvariableop/
+savev2_adam_v_11_kernel_read_readvariableop-
)savev2_adam_m_11_bias_read_readvariableop-
)savev2_adam_v_11_bias_read_readvariableop.
*savev2_adam_m_3_kernel_read_readvariableop.
*savev2_adam_v_3_kernel_read_readvariableop,
(savev2_adam_m_3_bias_read_readvariableop,
(savev2_adam_v_3_bias_read_readvariableop.
*savev2_adam_m_6_kernel_read_readvariableop.
*savev2_adam_v_6_kernel_read_readvariableop,
(savev2_adam_m_6_bias_read_readvariableop,
(savev2_adam_v_6_bias_read_readvariableop.
*savev2_adam_m_9_kernel_read_readvariableop.
*savev2_adam_v_9_kernel_read_readvariableop,
(savev2_adam_m_9_bias_read_readvariableop,
(savev2_adam_v_9_bias_read_readvariableop/
+savev2_adam_m_12_kernel_read_readvariableop/
+savev2_adam_v_12_kernel_read_readvariableop-
)savev2_adam_m_12_bias_read_readvariableop-
)savev2_adam_v_12_bias_read_readvariableop.
*savev2_adam_m_t_kernel_read_readvariableop.
*savev2_adam_v_t_kernel_read_readvariableop,
(savev2_adam_m_t_bias_read_readvariableop,
(savev2_adam_v_t_bias_read_readvariableop.
*savev2_adam_m_p_kernel_read_readvariableop.
*savev2_adam_v_p_kernel_read_readvariableop,
(savev2_adam_m_p_bias_read_readvariableop,
(savev2_adam_v_p_bias_read_readvariableop0
,savev2_adam_m_phi_kernel_read_readvariableop0
,savev2_adam_v_phi_kernel_read_readvariableop.
*savev2_adam_m_phi_bias_read_readvariableop.
*savev2_adam_v_phi_bias_read_readvariableop.
*savev2_adam_m_c_kernel_read_readvariableop.
*savev2_adam_v_c_kernel_read_readvariableop,
(savev2_adam_m_c_bias_read_readvariableop,
(savev2_adam_v_c_bias_read_readvariableop&
"savev2_total_4_read_readvariableop&
"savev2_count_4_read_readvariableop&
"savev2_total_3_read_readvariableop&
"savev2_count_3_read_readvariableop&
"savev2_total_2_read_readvariableop&
"savev2_count_2_read_readvariableop&
"savev2_total_1_read_readvariableop&
"savev2_count_1_read_readvariableop$
 savev2_total_read_readvariableop$
 savev2_count_read_readvariableop
savev2_const_2

identity_1��MergeV2Checkpointsw
StaticRegexFullMatchStaticRegexFullMatchfile_prefix"/device:CPU:**
_output_shapes
: *
pattern
^s3://.*Z
ConstConst"/device:CPU:**
_output_shapes
: *
dtype0*
valueB B.parta
Const_1Const"/device:CPU:**
_output_shapes
: *
dtype0*
valueB B
_temp/part�
SelectSelectStaticRegexFullMatch:output:0Const:output:0Const_1:output:0"/device:CPU:**
T0*
_output_shapes
: f

StringJoin
StringJoinfile_prefixSelect:output:0"/device:CPU:**
N*
_output_shapes
: L

num_shardsConst*
_output_shapes
: *
dtype0*
value	B :f
ShardedFilename/shardConst"/device:CPU:0*
_output_shapes
: *
dtype0*
value	B : �
ShardedFilenameShardedFilenameStringJoin:output:0ShardedFilename/shard:output:0num_shards:output:0"/device:CPU:0*
_output_shapes
: �6
SaveV2/tensor_namesConst"/device:CPU:0*
_output_shapes	
:�*
dtype0*�6
value�5B�5�B4layer_with_weights-0/mean/.ATTRIBUTES/VARIABLE_VALUEB8layer_with_weights-0/variance/.ATTRIBUTES/VARIABLE_VALUEB5layer_with_weights-0/count/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-3/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-3/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-4/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-4/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-5/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-5/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-6/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-6/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-7/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-7/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-8/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-8/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-9/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-9/bias/.ATTRIBUTES/VARIABLE_VALUEB7layer_with_weights-10/kernel/.ATTRIBUTES/VARIABLE_VALUEB5layer_with_weights-10/bias/.ATTRIBUTES/VARIABLE_VALUEB7layer_with_weights-11/kernel/.ATTRIBUTES/VARIABLE_VALUEB5layer_with_weights-11/bias/.ATTRIBUTES/VARIABLE_VALUEB7layer_with_weights-12/kernel/.ATTRIBUTES/VARIABLE_VALUEB5layer_with_weights-12/bias/.ATTRIBUTES/VARIABLE_VALUEB7layer_with_weights-13/kernel/.ATTRIBUTES/VARIABLE_VALUEB5layer_with_weights-13/bias/.ATTRIBUTES/VARIABLE_VALUEB7layer_with_weights-14/kernel/.ATTRIBUTES/VARIABLE_VALUEB5layer_with_weights-14/bias/.ATTRIBUTES/VARIABLE_VALUEB7layer_with_weights-15/kernel/.ATTRIBUTES/VARIABLE_VALUEB5layer_with_weights-15/bias/.ATTRIBUTES/VARIABLE_VALUEB7layer_with_weights-16/kernel/.ATTRIBUTES/VARIABLE_VALUEB5layer_with_weights-16/bias/.ATTRIBUTES/VARIABLE_VALUEB7layer_with_weights-17/kernel/.ATTRIBUTES/VARIABLE_VALUEB5layer_with_weights-17/bias/.ATTRIBUTES/VARIABLE_VALUEB7layer_with_weights-18/kernel/.ATTRIBUTES/VARIABLE_VALUEB5layer_with_weights-18/bias/.ATTRIBUTES/VARIABLE_VALUEB7layer_with_weights-19/kernel/.ATTRIBUTES/VARIABLE_VALUEB5layer_with_weights-19/bias/.ATTRIBUTES/VARIABLE_VALUEB0optimizer/_iterations/.ATTRIBUTES/VARIABLE_VALUEB3optimizer/_learning_rate/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/1/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/2/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/3/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/4/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/5/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/6/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/7/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/8/.ATTRIBUTES/VARIABLE_VALUEB1optimizer/_variables/9/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/10/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/11/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/12/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/13/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/14/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/15/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/16/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/17/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/18/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/19/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/20/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/21/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/22/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/23/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/24/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/25/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/26/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/27/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/28/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/29/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/30/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/31/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/32/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/33/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/34/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/35/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/36/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/37/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/38/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/39/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/40/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/41/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/42/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/43/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/44/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/45/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/46/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/47/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/48/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/49/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/50/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/51/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/52/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/53/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/54/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/55/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/56/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/57/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/58/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/59/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/60/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/61/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/62/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/63/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/64/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/65/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/66/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/67/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/68/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/69/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/70/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/71/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/72/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/73/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/74/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/75/.ATTRIBUTES/VARIABLE_VALUEB2optimizer/_variables/76/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/0/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/1/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/2/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/2/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/3/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/3/count/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/4/total/.ATTRIBUTES/VARIABLE_VALUEB4keras_api/metrics/4/count/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH�
SaveV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes	
:�*
dtype0*�
value�B��B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B B �+
SaveV2SaveV2ShardedFilename:filename:0SaveV2/tensor_names:output:0 SaveV2/shape_and_slices:output:0savev2_mean_read_readvariableop#savev2_variance_read_readvariableop"savev2_count_5_read_readvariableop$savev2_l2_kernel_read_readvariableop"savev2_l2_bias_read_readvariableop$savev2_l3_kernel_read_readvariableop"savev2_l3_bias_read_readvariableop#savev2_x_kernel_read_readvariableop!savev2_x_bias_read_readvariableop#savev2_1_kernel_read_readvariableop!savev2_1_bias_read_readvariableop#savev2_4_kernel_read_readvariableop!savev2_4_bias_read_readvariableop#savev2_7_kernel_read_readvariableop!savev2_7_bias_read_readvariableop$savev2_10_kernel_read_readvariableop"savev2_10_bias_read_readvariableop#savev2_2_kernel_read_readvariableop!savev2_2_bias_read_readvariableop#savev2_5_kernel_read_readvariableop!savev2_5_bias_read_readvariableop#savev2_8_kernel_read_readvariableop!savev2_8_bias_read_readvariableop$savev2_11_kernel_read_readvariableop"savev2_11_bias_read_readvariableop#savev2_3_kernel_read_readvariableop!savev2_3_bias_read_readvariableop#savev2_6_kernel_read_readvariableop!savev2_6_bias_read_readvariableop#savev2_9_kernel_read_readvariableop!savev2_9_bias_read_readvariableop$savev2_12_kernel_read_readvariableop"savev2_12_bias_read_readvariableop#savev2_t_kernel_read_readvariableop!savev2_t_bias_read_readvariableop#savev2_p_kernel_read_readvariableop!savev2_p_bias_read_readvariableop%savev2_phi_kernel_read_readvariableop#savev2_phi_bias_read_readvariableop#savev2_c_kernel_read_readvariableop!savev2_c_bias_read_readvariableop$savev2_iteration_read_readvariableop(savev2_learning_rate_read_readvariableop+savev2_adam_m_l2_kernel_read_readvariableop+savev2_adam_v_l2_kernel_read_readvariableop)savev2_adam_m_l2_bias_read_readvariableop)savev2_adam_v_l2_bias_read_readvariableop+savev2_adam_m_l3_kernel_read_readvariableop+savev2_adam_v_l3_kernel_read_readvariableop)savev2_adam_m_l3_bias_read_readvariableop)savev2_adam_v_l3_bias_read_readvariableop*savev2_adam_m_x_kernel_read_readvariableop*savev2_adam_v_x_kernel_read_readvariableop(savev2_adam_m_x_bias_read_readvariableop(savev2_adam_v_x_bias_read_readvariableop*savev2_adam_m_1_kernel_read_readvariableop*savev2_adam_v_1_kernel_read_readvariableop(savev2_adam_m_1_bias_read_readvariableop(savev2_adam_v_1_bias_read_readvariableop*savev2_adam_m_4_kernel_read_readvariableop*savev2_adam_v_4_kernel_read_readvariableop(savev2_adam_m_4_bias_read_readvariableop(savev2_adam_v_4_bias_read_readvariableop*savev2_adam_m_7_kernel_read_readvariableop*savev2_adam_v_7_kernel_read_readvariableop(savev2_adam_m_7_bias_read_readvariableop(savev2_adam_v_7_bias_read_readvariableop+savev2_adam_m_10_kernel_read_readvariableop+savev2_adam_v_10_kernel_read_readvariableop)savev2_adam_m_10_bias_read_readvariableop)savev2_adam_v_10_bias_read_readvariableop*savev2_adam_m_2_kernel_read_readvariableop*savev2_adam_v_2_kernel_read_readvariableop(savev2_adam_m_2_bias_read_readvariableop(savev2_adam_v_2_bias_read_readvariableop*savev2_adam_m_5_kernel_read_readvariableop*savev2_adam_v_5_kernel_read_readvariableop(savev2_adam_m_5_bias_read_readvariableop(savev2_adam_v_5_bias_read_readvariableop*savev2_adam_m_8_kernel_read_readvariableop*savev2_adam_v_8_kernel_read_readvariableop(savev2_adam_m_8_bias_read_readvariableop(savev2_adam_v_8_bias_read_readvariableop+savev2_adam_m_11_kernel_read_readvariableop+savev2_adam_v_11_kernel_read_readvariableop)savev2_adam_m_11_bias_read_readvariableop)savev2_adam_v_11_bias_read_readvariableop*savev2_adam_m_3_kernel_read_readvariableop*savev2_adam_v_3_kernel_read_readvariableop(savev2_adam_m_3_bias_read_readvariableop(savev2_adam_v_3_bias_read_readvariableop*savev2_adam_m_6_kernel_read_readvariableop*savev2_adam_v_6_kernel_read_readvariableop(savev2_adam_m_6_bias_read_readvariableop(savev2_adam_v_6_bias_read_readvariableop*savev2_adam_m_9_kernel_read_readvariableop*savev2_adam_v_9_kernel_read_readvariableop(savev2_adam_m_9_bias_read_readvariableop(savev2_adam_v_9_bias_read_readvariableop+savev2_adam_m_12_kernel_read_readvariableop+savev2_adam_v_12_kernel_read_readvariableop)savev2_adam_m_12_bias_read_readvariableop)savev2_adam_v_12_bias_read_readvariableop*savev2_adam_m_t_kernel_read_readvariableop*savev2_adam_v_t_kernel_read_readvariableop(savev2_adam_m_t_bias_read_readvariableop(savev2_adam_v_t_bias_read_readvariableop*savev2_adam_m_p_kernel_read_readvariableop*savev2_adam_v_p_kernel_read_readvariableop(savev2_adam_m_p_bias_read_readvariableop(savev2_adam_v_p_bias_read_readvariableop,savev2_adam_m_phi_kernel_read_readvariableop,savev2_adam_v_phi_kernel_read_readvariableop*savev2_adam_m_phi_bias_read_readvariableop*savev2_adam_v_phi_bias_read_readvariableop*savev2_adam_m_c_kernel_read_readvariableop*savev2_adam_v_c_kernel_read_readvariableop(savev2_adam_m_c_bias_read_readvariableop(savev2_adam_v_c_bias_read_readvariableop"savev2_total_4_read_readvariableop"savev2_count_4_read_readvariableop"savev2_total_3_read_readvariableop"savev2_count_3_read_readvariableop"savev2_total_2_read_readvariableop"savev2_count_2_read_readvariableop"savev2_total_1_read_readvariableop"savev2_count_1_read_readvariableop savev2_total_read_readvariableop savev2_count_read_readvariableopsavev2_const_2"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 *�
dtypes�
�2�		�
&MergeV2Checkpoints/checkpoint_prefixesPackShardedFilename:filename:0^SaveV2"/device:CPU:0*
N*
T0*
_output_shapes
:�
MergeV2CheckpointsMergeV2Checkpoints/MergeV2Checkpoints/checkpoint_prefixes:output:0file_prefix"/device:CPU:0*&
 _has_manual_control_dependencies(*
_output_shapes
 f
IdentityIdentityfile_prefix^MergeV2Checkpoints"/device:CPU:0*
T0*
_output_shapes
: Q

Identity_1IdentityIdentity:output:0^NoOp*
T0*
_output_shapes
: [
NoOpNoOp^MergeV2Checkpoints*"
_acd_function_control_output(*
_output_shapes
 "!

identity_1Identity_1:output:0*�
_input_shapes�
�: ::: :@:@:@@:@:@@:@:@@:@:@@:@:@@:@:@@:@:@@:@:@@:@:@@:@:@@:@:@ : :@ : :@ : :@ : : :: :: :: :: : :@:@:@:@:@@:@@:@:@:@@:@@:@:@:@@:@@:@:@:@@:@@:@:@:@@:@@:@:@:@@:@@:@:@:@@:@@:@:@:@@:@@:@:@:@@:@@:@:@:@@:@@:@:@:@ :@ : : :@ :@ : : :@ :@ : : :@ :@ : : : : ::: : ::: : ::: : ::: : : : : : : : : : : 2(
MergeV2CheckpointsMergeV2Checkpoints:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix: 

_output_shapes
:: 

_output_shapes
::

_output_shapes
: :$ 

_output_shapes

:@: 

_output_shapes
:@:$ 

_output_shapes

:@@: 

_output_shapes
:@:$ 

_output_shapes

:@@: 	

_output_shapes
:@:$
 

_output_shapes

:@@: 

_output_shapes
:@:$ 

_output_shapes

:@@: 

_output_shapes
:@:$ 

_output_shapes

:@@: 

_output_shapes
:@:$ 

_output_shapes

:@@: 

_output_shapes
:@:$ 

_output_shapes

:@@: 

_output_shapes
:@:$ 

_output_shapes

:@@: 

_output_shapes
:@:$ 

_output_shapes

:@@: 

_output_shapes
:@:$ 

_output_shapes

:@@: 

_output_shapes
:@:$ 

_output_shapes

:@ : 

_output_shapes
: :$ 

_output_shapes

:@ : 

_output_shapes
: :$ 

_output_shapes

:@ : 

_output_shapes
: :$  

_output_shapes

:@ : !

_output_shapes
: :$" 

_output_shapes

: : #

_output_shapes
::$$ 

_output_shapes

: : %

_output_shapes
::$& 

_output_shapes

: : '

_output_shapes
::$( 

_output_shapes

: : )

_output_shapes
::*

_output_shapes
: :+

_output_shapes
: :$, 

_output_shapes

:@:$- 

_output_shapes

:@: .

_output_shapes
:@: /

_output_shapes
:@:$0 

_output_shapes

:@@:$1 

_output_shapes

:@@: 2

_output_shapes
:@: 3

_output_shapes
:@:$4 

_output_shapes

:@@:$5 

_output_shapes

:@@: 6

_output_shapes
:@: 7

_output_shapes
:@:$8 

_output_shapes

:@@:$9 

_output_shapes

:@@: :

_output_shapes
:@: ;

_output_shapes
:@:$< 

_output_shapes

:@@:$= 

_output_shapes

:@@: >

_output_shapes
:@: ?

_output_shapes
:@:$@ 

_output_shapes

:@@:$A 

_output_shapes

:@@: B

_output_shapes
:@: C

_output_shapes
:@:$D 

_output_shapes

:@@:$E 

_output_shapes

:@@: F

_output_shapes
:@: G

_output_shapes
:@:$H 

_output_shapes

:@@:$I 

_output_shapes

:@@: J

_output_shapes
:@: K

_output_shapes
:@:$L 

_output_shapes

:@@:$M 

_output_shapes

:@@: N

_output_shapes
:@: O

_output_shapes
:@:$P 

_output_shapes

:@@:$Q 

_output_shapes

:@@: R

_output_shapes
:@: S

_output_shapes
:@:$T 

_output_shapes

:@@:$U 

_output_shapes

:@@: V

_output_shapes
:@: W

_output_shapes
:@:$X 

_output_shapes

:@ :$Y 

_output_shapes

:@ : Z

_output_shapes
: : [

_output_shapes
: :$\ 

_output_shapes

:@ :$] 

_output_shapes

:@ : ^

_output_shapes
: : _

_output_shapes
: :$` 

_output_shapes

:@ :$a 

_output_shapes

:@ : b

_output_shapes
: : c

_output_shapes
: :$d 

_output_shapes

:@ :$e 

_output_shapes

:@ : f

_output_shapes
: : g

_output_shapes
: :$h 

_output_shapes

: :$i 

_output_shapes

: : j

_output_shapes
:: k

_output_shapes
::$l 

_output_shapes

: :$m 

_output_shapes

: : n

_output_shapes
:: o

_output_shapes
::$p 

_output_shapes

: :$q 

_output_shapes

: : r

_output_shapes
:: s

_output_shapes
::$t 

_output_shapes

: :$u 

_output_shapes

: : v

_output_shapes
:: w

_output_shapes
::x

_output_shapes
: :y

_output_shapes
: :z

_output_shapes
: :{

_output_shapes
: :|

_output_shapes
: :}

_output_shapes
: :~

_output_shapes
: :

_output_shapes
: :�

_output_shapes
: :�

_output_shapes
: :�

_output_shapes
: 
�
�
$__inference_l2_layer_call_fn_2218709

inputs
unknown:@
	unknown_0:@
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *H
fCRA
?__inference_l2_layer_call_and_return_conditional_losses_2216889o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������@`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�	
'__inference_model_layer_call_fn_2218315

inputs
unknown
	unknown_0
	unknown_1:@
	unknown_2:@
	unknown_3:@@
	unknown_4:@
	unknown_5:@@
	unknown_6:@
	unknown_7:@@
	unknown_8:@
	unknown_9:@@

unknown_10:@

unknown_11:@@

unknown_12:@

unknown_13:@@

unknown_14:@

unknown_15:@@

unknown_16:@

unknown_17:@@

unknown_18:@

unknown_19:@@

unknown_20:@

unknown_21:@@

unknown_22:@

unknown_23:@ 

unknown_24: 

unknown_25:@ 

unknown_26: 

unknown_27:@ 

unknown_28: 

unknown_29:@ 

unknown_30: 

unknown_31: 

unknown_32:

unknown_33: 

unknown_34:

unknown_35: 

unknown_36:

unknown_37: 

unknown_38:
identity

identity_1

identity_2

identity_3��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17
unknown_18
unknown_19
unknown_20
unknown_21
unknown_22
unknown_23
unknown_24
unknown_25
unknown_26
unknown_27
unknown_28
unknown_29
unknown_30
unknown_31
unknown_32
unknown_33
unknown_34
unknown_35
unknown_36
unknown_37
unknown_38*4
Tin-
+2)*
Tout
2*
_collective_manager_ids
 *`
_output_shapesN
L:���������:���������:���������:���������*H
_read_only_resource_inputs*
(&	
 !"#$%&'(*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_model_layer_call_and_return_conditional_losses_2217205o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������q

Identity_1Identity StatefulPartitionedCall:output:1^NoOp*
T0*'
_output_shapes
:���������q

Identity_2Identity StatefulPartitionedCall:output:2^NoOp*
T0*'
_output_shapes
:���������q

Identity_3Identity StatefulPartitionedCall:output:3^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0"!

identity_3Identity_3:output:0*(
_construction_contextkEagerRuntime*�
_input_shapesu
s:���������::: : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs:$ 

_output_shapes

::$ 

_output_shapes

:
�

�
?__inference_10_layer_call_and_return_conditional_losses_2216940

inputs0
matmul_readvariableop_resource:@@-
biasadd_readvariableop_resource:@
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@@*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������@a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:���������@w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�
�	
'__inference_model_layer_call_fn_2218406

inputs
unknown
	unknown_0
	unknown_1:@
	unknown_2:@
	unknown_3:@@
	unknown_4:@
	unknown_5:@@
	unknown_6:@
	unknown_7:@@
	unknown_8:@
	unknown_9:@@

unknown_10:@

unknown_11:@@

unknown_12:@

unknown_13:@@

unknown_14:@

unknown_15:@@

unknown_16:@

unknown_17:@@

unknown_18:@

unknown_19:@@

unknown_20:@

unknown_21:@@

unknown_22:@

unknown_23:@ 

unknown_24: 

unknown_25:@ 

unknown_26: 

unknown_27:@ 

unknown_28: 

unknown_29:@ 

unknown_30: 

unknown_31: 

unknown_32:

unknown_33: 

unknown_34:

unknown_35: 

unknown_36:

unknown_37: 

unknown_38:
identity

identity_1

identity_2

identity_3��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17
unknown_18
unknown_19
unknown_20
unknown_21
unknown_22
unknown_23
unknown_24
unknown_25
unknown_26
unknown_27
unknown_28
unknown_29
unknown_30
unknown_31
unknown_32
unknown_33
unknown_34
unknown_35
unknown_36
unknown_37
unknown_38*4
Tin-
+2)*
Tout
2*
_collective_manager_ids
 *`
_output_shapesN
L:���������:���������:���������:���������*H
_read_only_resource_inputs*
(&	
 !"#$%&'(*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_model_layer_call_and_return_conditional_losses_2217686o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������q

Identity_1Identity StatefulPartitionedCall:output:1^NoOp*
T0*'
_output_shapes
:���������q

Identity_2Identity StatefulPartitionedCall:output:2^NoOp*
T0*'
_output_shapes
:���������q

Identity_3Identity StatefulPartitionedCall:output:3^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0"!

identity_3Identity_3:output:0*(
_construction_contextkEagerRuntime*�
_input_shapesu
s:���������::: : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs:$ 

_output_shapes

::$ 

_output_shapes

:
�^
�
B__inference_model_layer_call_and_return_conditional_losses_2218084	
input
normalization_sub_y
normalization_sqrt_x

l2_2217985:@

l2_2217987:@

l3_2217990:@@

l3_2217992:@
	x_2217995:@@
	x_2217997:@
unknown:@@
	unknown_0:@
	unknown_1:@@
	unknown_2:@
	unknown_3:@@
	unknown_4:@
	unknown_5:@@
	unknown_6:@
	unknown_7:@@
	unknown_8:@
	unknown_9:@@

unknown_10:@

unknown_11:@@

unknown_12:@

unknown_13:@@

unknown_14:@

unknown_15:@ 

unknown_16: 

unknown_17:@ 

unknown_18: 

unknown_19:@ 

unknown_20: 

unknown_21:@ 

unknown_22: 
	c_2218060: 
	c_2218062:
phi_2218065: 
phi_2218067:
	p_2218070: 
	p_2218072:
	t_2218075: 
	t_2218077:
identity

identity_1

identity_2

identity_3��1/StatefulPartitionedCall�10/StatefulPartitionedCall�11/StatefulPartitionedCall�12/StatefulPartitionedCall�2/StatefulPartitionedCall�3/StatefulPartitionedCall�4/StatefulPartitionedCall�5/StatefulPartitionedCall�6/StatefulPartitionedCall�7/StatefulPartitionedCall�8/StatefulPartitionedCall�9/StatefulPartitionedCall�P/StatefulPartitionedCall�T/StatefulPartitionedCall�c/StatefulPartitionedCall�l2/StatefulPartitionedCall�l3/StatefulPartitionedCall�phi/StatefulPartitionedCall�x/StatefulPartitionedCallf
normalization/subSubinputnormalization_sub_y*
T0*'
_output_shapes
:���������Y
normalization/SqrtSqrtnormalization_sqrt_x*
T0*
_output_shapes

:\
normalization/Maximum/yConst*
_output_shapes
: *
dtype0*
valueB
 *���3�
normalization/MaximumMaximumnormalization/Sqrt:y:0 normalization/Maximum/y:output:0*
T0*
_output_shapes

:�
normalization/truedivRealDivnormalization/sub:z:0normalization/Maximum:z:0*
T0*'
_output_shapes
:����������
l2/StatefulPartitionedCallStatefulPartitionedCallnormalization/truediv:z:0
l2_2217985
l2_2217987*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *H
fCRA
?__inference_l2_layer_call_and_return_conditional_losses_2216889�
l3/StatefulPartitionedCallStatefulPartitionedCall#l2/StatefulPartitionedCall:output:0
l3_2217990
l3_2217992*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *H
fCRA
?__inference_l3_layer_call_and_return_conditional_losses_2216906�
x/StatefulPartitionedCallStatefulPartitionedCall#l3/StatefulPartitionedCall:output:0	x_2217995	x_2217997*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_x_layer_call_and_return_conditional_losses_2216923�
10/StatefulPartitionedCallStatefulPartitionedCall"x/StatefulPartitionedCall:output:0unknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *H
fCRA
?__inference_10_layer_call_and_return_conditional_losses_2216940�
7/StatefulPartitionedCallStatefulPartitionedCall"x/StatefulPartitionedCall:output:0	unknown_1	unknown_2*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_7_layer_call_and_return_conditional_losses_2216957�
4/StatefulPartitionedCallStatefulPartitionedCall"x/StatefulPartitionedCall:output:0	unknown_3	unknown_4*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_4_layer_call_and_return_conditional_losses_2216974�
1/StatefulPartitionedCallStatefulPartitionedCall"x/StatefulPartitionedCall:output:0	unknown_5	unknown_6*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_1_layer_call_and_return_conditional_losses_2216991�
11/StatefulPartitionedCallStatefulPartitionedCall#10/StatefulPartitionedCall:output:0	unknown_7	unknown_8*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *H
fCRA
?__inference_11_layer_call_and_return_conditional_losses_2217008�
8/StatefulPartitionedCallStatefulPartitionedCall"7/StatefulPartitionedCall:output:0	unknown_9
unknown_10*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_8_layer_call_and_return_conditional_losses_2217025�
5/StatefulPartitionedCallStatefulPartitionedCall"4/StatefulPartitionedCall:output:0
unknown_11
unknown_12*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_5_layer_call_and_return_conditional_losses_2217042�
2/StatefulPartitionedCallStatefulPartitionedCall"1/StatefulPartitionedCall:output:0
unknown_13
unknown_14*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_2_layer_call_and_return_conditional_losses_2217059�
12/StatefulPartitionedCallStatefulPartitionedCall#11/StatefulPartitionedCall:output:0
unknown_15
unknown_16*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:��������� *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *H
fCRA
?__inference_12_layer_call_and_return_conditional_losses_2217076�
9/StatefulPartitionedCallStatefulPartitionedCall"8/StatefulPartitionedCall:output:0
unknown_17
unknown_18*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:��������� *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_9_layer_call_and_return_conditional_losses_2217093�
6/StatefulPartitionedCallStatefulPartitionedCall"5/StatefulPartitionedCall:output:0
unknown_19
unknown_20*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:��������� *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_6_layer_call_and_return_conditional_losses_2217110�
3/StatefulPartitionedCallStatefulPartitionedCall"2/StatefulPartitionedCall:output:0
unknown_21
unknown_22*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:��������� *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_3_layer_call_and_return_conditional_losses_2217127�
c/StatefulPartitionedCallStatefulPartitionedCall#12/StatefulPartitionedCall:output:0	c_2218060	c_2218062*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_c_layer_call_and_return_conditional_losses_2217144�
phi/StatefulPartitionedCallStatefulPartitionedCall"9/StatefulPartitionedCall:output:0phi_2218065phi_2218067*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *I
fDRB
@__inference_phi_layer_call_and_return_conditional_losses_2217161�
P/StatefulPartitionedCallStatefulPartitionedCall"6/StatefulPartitionedCall:output:0	p_2218070	p_2218072*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_P_layer_call_and_return_conditional_losses_2217178�
T/StatefulPartitionedCallStatefulPartitionedCall"3/StatefulPartitionedCall:output:0	t_2218075	t_2218077*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_T_layer_call_and_return_conditional_losses_2217195q
IdentityIdentity"T/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������s

Identity_1Identity"P/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������u

Identity_2Identity$phi/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������s

Identity_3Identity"c/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:����������
NoOpNoOp^1/StatefulPartitionedCall^10/StatefulPartitionedCall^11/StatefulPartitionedCall^12/StatefulPartitionedCall^2/StatefulPartitionedCall^3/StatefulPartitionedCall^4/StatefulPartitionedCall^5/StatefulPartitionedCall^6/StatefulPartitionedCall^7/StatefulPartitionedCall^8/StatefulPartitionedCall^9/StatefulPartitionedCall^P/StatefulPartitionedCall^T/StatefulPartitionedCall^c/StatefulPartitionedCall^l2/StatefulPartitionedCall^l3/StatefulPartitionedCall^phi/StatefulPartitionedCall^x/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0"!

identity_3Identity_3:output:0*(
_construction_contextkEagerRuntime*�
_input_shapesu
s:���������::: : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : 26
1/StatefulPartitionedCall1/StatefulPartitionedCall28
10/StatefulPartitionedCall10/StatefulPartitionedCall28
11/StatefulPartitionedCall11/StatefulPartitionedCall28
12/StatefulPartitionedCall12/StatefulPartitionedCall26
2/StatefulPartitionedCall2/StatefulPartitionedCall26
3/StatefulPartitionedCall3/StatefulPartitionedCall26
4/StatefulPartitionedCall4/StatefulPartitionedCall26
5/StatefulPartitionedCall5/StatefulPartitionedCall26
6/StatefulPartitionedCall6/StatefulPartitionedCall26
7/StatefulPartitionedCall7/StatefulPartitionedCall26
8/StatefulPartitionedCall8/StatefulPartitionedCall26
9/StatefulPartitionedCall9/StatefulPartitionedCall26
P/StatefulPartitionedCallP/StatefulPartitionedCall26
T/StatefulPartitionedCallT/StatefulPartitionedCall26
c/StatefulPartitionedCallc/StatefulPartitionedCall28
l2/StatefulPartitionedCalll2/StatefulPartitionedCall28
l3/StatefulPartitionedCalll3/StatefulPartitionedCall2:
phi/StatefulPartitionedCallphi/StatefulPartitionedCall26
x/StatefulPartitionedCallx/StatefulPartitionedCall:N J
'
_output_shapes
:���������

_user_specified_nameinput:$ 

_output_shapes

::$ 

_output_shapes

:
�

�
>__inference_6_layer_call_and_return_conditional_losses_2218960

inputs0
matmul_readvariableop_resource:@ -
biasadd_readvariableop_resource: 
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@ *
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:��������� a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:��������� w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�
�
#__inference_3_layer_call_fn_2218929

inputs
unknown:@ 
	unknown_0: 
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:��������� *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_3_layer_call_and_return_conditional_losses_2217127o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:��������� `
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�
�
#__inference_x_layer_call_fn_2218749

inputs
unknown:@@
	unknown_0:@
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_x_layer_call_and_return_conditional_losses_2216923o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������@`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�
�
#__inference_2_layer_call_fn_2218849

inputs
unknown:@@
	unknown_0:@
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_2_layer_call_and_return_conditional_losses_2217059o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������@`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�

�
?__inference_l2_layer_call_and_return_conditional_losses_2218720

inputs0
matmul_readvariableop_resource:@-
biasadd_readvariableop_resource:@
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������@a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:���������@w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�

�
?__inference_12_layer_call_and_return_conditional_losses_2219000

inputs0
matmul_readvariableop_resource:@ -
biasadd_readvariableop_resource: 
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@ *
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:��������� a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:��������� w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�
�	
'__inference_model_layer_call_fn_2217866	
input
unknown
	unknown_0
	unknown_1:@
	unknown_2:@
	unknown_3:@@
	unknown_4:@
	unknown_5:@@
	unknown_6:@
	unknown_7:@@
	unknown_8:@
	unknown_9:@@

unknown_10:@

unknown_11:@@

unknown_12:@

unknown_13:@@

unknown_14:@

unknown_15:@@

unknown_16:@

unknown_17:@@

unknown_18:@

unknown_19:@@

unknown_20:@

unknown_21:@@

unknown_22:@

unknown_23:@ 

unknown_24: 

unknown_25:@ 

unknown_26: 

unknown_27:@ 

unknown_28: 

unknown_29:@ 

unknown_30: 

unknown_31: 

unknown_32:

unknown_33: 

unknown_34:

unknown_35: 

unknown_36:

unknown_37: 

unknown_38:
identity

identity_1

identity_2

identity_3��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17
unknown_18
unknown_19
unknown_20
unknown_21
unknown_22
unknown_23
unknown_24
unknown_25
unknown_26
unknown_27
unknown_28
unknown_29
unknown_30
unknown_31
unknown_32
unknown_33
unknown_34
unknown_35
unknown_36
unknown_37
unknown_38*4
Tin-
+2)*
Tout
2*
_collective_manager_ids
 *`
_output_shapesN
L:���������:���������:���������:���������*H
_read_only_resource_inputs*
(&	
 !"#$%&'(*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_model_layer_call_and_return_conditional_losses_2217686o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������q

Identity_1Identity StatefulPartitionedCall:output:1^NoOp*
T0*'
_output_shapes
:���������q

Identity_2Identity StatefulPartitionedCall:output:2^NoOp*
T0*'
_output_shapes
:���������q

Identity_3Identity StatefulPartitionedCall:output:3^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0"!

identity_3Identity_3:output:0*(
_construction_contextkEagerRuntime*�
_input_shapesu
s:���������::: : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:N J
'
_output_shapes
:���������

_user_specified_nameinput:$ 

_output_shapes

::$ 

_output_shapes

:
�

�
?__inference_12_layer_call_and_return_conditional_losses_2217076

inputs0
matmul_readvariableop_resource:@ -
biasadd_readvariableop_resource: 
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@ *
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:��������� a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:��������� w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�
�
#__inference_4_layer_call_fn_2218789

inputs
unknown:@@
	unknown_0:@
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_4_layer_call_and_return_conditional_losses_2216974o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������@`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�

�
>__inference_c_layer_call_and_return_conditional_losses_2219080

inputs0
matmul_readvariableop_resource: -
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

: *
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:���������w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:��������� : : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:��������� 
 
_user_specified_nameinputs
�

�
?__inference_l2_layer_call_and_return_conditional_losses_2216889

inputs0
matmul_readvariableop_resource:@-
biasadd_readvariableop_resource:@
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������@a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:���������@w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������
 
_user_specified_nameinputs
�
�
#__inference_c_layer_call_fn_2219069

inputs
unknown: 
	unknown_0:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_c_layer_call_and_return_conditional_losses_2217144o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:��������� : : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:��������� 
 
_user_specified_nameinputs
�

�
>__inference_1_layer_call_and_return_conditional_losses_2218780

inputs0
matmul_readvariableop_resource:@@-
biasadd_readvariableop_resource:@
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@@*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������@a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:���������@w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�
�
%__inference_phi_layer_call_fn_2219049

inputs
unknown: 
	unknown_0:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *I
fDRB
@__inference_phi_layer_call_and_return_conditional_losses_2217161o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:��������� : : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:��������� 
 
_user_specified_nameinputs
�

�
>__inference_5_layer_call_and_return_conditional_losses_2218880

inputs0
matmul_readvariableop_resource:@@-
biasadd_readvariableop_resource:@
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@@*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������@a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:���������@w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�
�
$__inference_12_layer_call_fn_2218989

inputs
unknown:@ 
	unknown_0: 
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:��������� *$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *H
fCRA
?__inference_12_layer_call_and_return_conditional_losses_2217076o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:��������� `
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�

�
>__inference_9_layer_call_and_return_conditional_losses_2217093

inputs0
matmul_readvariableop_resource:@ -
biasadd_readvariableop_resource: 
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@ *
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
: *
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:��������� P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:��������� a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:��������� w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�

�
?__inference_l3_layer_call_and_return_conditional_losses_2218740

inputs0
matmul_readvariableop_resource:@@-
biasadd_readvariableop_resource:@
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@@*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������@a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:���������@w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�

�
?__inference_11_layer_call_and_return_conditional_losses_2218920

inputs0
matmul_readvariableop_resource:@@-
biasadd_readvariableop_resource:@
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@@*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������@a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:���������@w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�
�
#__inference_1_layer_call_fn_2218769

inputs
unknown:@@
	unknown_0:@
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_1_layer_call_and_return_conditional_losses_2216991o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������@`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�
�	
'__inference_model_layer_call_fn_2217294	
input
unknown
	unknown_0
	unknown_1:@
	unknown_2:@
	unknown_3:@@
	unknown_4:@
	unknown_5:@@
	unknown_6:@
	unknown_7:@@
	unknown_8:@
	unknown_9:@@

unknown_10:@

unknown_11:@@

unknown_12:@

unknown_13:@@

unknown_14:@

unknown_15:@@

unknown_16:@

unknown_17:@@

unknown_18:@

unknown_19:@@

unknown_20:@

unknown_21:@@

unknown_22:@

unknown_23:@ 

unknown_24: 

unknown_25:@ 

unknown_26: 

unknown_27:@ 

unknown_28: 

unknown_29:@ 

unknown_30: 

unknown_31: 

unknown_32:

unknown_33: 

unknown_34:

unknown_35: 

unknown_36:

unknown_37: 

unknown_38:
identity

identity_1

identity_2

identity_3��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17
unknown_18
unknown_19
unknown_20
unknown_21
unknown_22
unknown_23
unknown_24
unknown_25
unknown_26
unknown_27
unknown_28
unknown_29
unknown_30
unknown_31
unknown_32
unknown_33
unknown_34
unknown_35
unknown_36
unknown_37
unknown_38*4
Tin-
+2)*
Tout
2*
_collective_manager_ids
 *`
_output_shapesN
L:���������:���������:���������:���������*H
_read_only_resource_inputs*
(&	
 !"#$%&'(*-
config_proto

CPU

GPU 2J 8� *K
fFRD
B__inference_model_layer_call_and_return_conditional_losses_2217205o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������q

Identity_1Identity StatefulPartitionedCall:output:1^NoOp*
T0*'
_output_shapes
:���������q

Identity_2Identity StatefulPartitionedCall:output:2^NoOp*
T0*'
_output_shapes
:���������q

Identity_3Identity StatefulPartitionedCall:output:3^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0"!

identity_1Identity_1:output:0"!

identity_2Identity_2:output:0"!

identity_3Identity_3:output:0*(
_construction_contextkEagerRuntime*�
_input_shapesu
s:���������::: : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:N J
'
_output_shapes
:���������

_user_specified_nameinput:$ 

_output_shapes

::$ 

_output_shapes

:
�
�
#__inference_5_layer_call_fn_2218869

inputs
unknown:@@
	unknown_0:@
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_5_layer_call_and_return_conditional_losses_2217042o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������@`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�
�
#__inference_P_layer_call_fn_2219029

inputs
unknown: 
	unknown_0:
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_P_layer_call_and_return_conditional_losses_2217178o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:��������� : : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:��������� 
 
_user_specified_nameinputs
�

�
>__inference_P_layer_call_and_return_conditional_losses_2217178

inputs0
matmul_readvariableop_resource: -
biasadd_readvariableop_resource:
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

: *
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:���������w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:��������� : : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:��������� 
 
_user_specified_nameinputs
�

�
>__inference_x_layer_call_and_return_conditional_losses_2218760

inputs0
matmul_readvariableop_resource:@@-
biasadd_readvariableop_resource:@
identity��BiasAdd/ReadVariableOp�MatMul/ReadVariableOpt
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes

:@@*
dtype0i
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@r
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype0v
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:���������@P
ReluReluBiasAdd:output:0*
T0*'
_output_shapes
:���������@a
IdentityIdentityRelu:activations:0^NoOp*
T0*'
_output_shapes
:���������@w
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs
�
�
#__inference_7_layer_call_fn_2218809

inputs
unknown:@@
	unknown_0:@
identity��StatefulPartitionedCall�
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:���������@*$
_read_only_resource_inputs
*-
config_proto

CPU

GPU 2J 8� *G
fBR@
>__inference_7_layer_call_and_return_conditional_losses_2216957o
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:���������@`
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 "
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:���������@: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:���������@
 
_user_specified_nameinputs"�
L
saver_filename:0StatefulPartitionedCall_1:0StatefulPartitionedCall_28"
saved_model_main_op

NoOp*>
__saved_model_init_op%#
__saved_model_init_op

NoOp*�
serving_default�
7
input.
serving_default_input:0���������5
P0
StatefulPartitionedCall:0���������5
T0
StatefulPartitionedCall:1���������5
c0
StatefulPartitionedCall:2���������7
phi0
StatefulPartitionedCall:3���������tensorflow/serving/predict:��
�
layer-0
layer_with_weights-0
layer-1
layer_with_weights-1
layer-2
layer_with_weights-2
layer-3
layer_with_weights-3
layer-4
layer_with_weights-4
layer-5
layer_with_weights-5
layer-6
layer_with_weights-6
layer-7
	layer_with_weights-7
	layer-8

layer_with_weights-8

layer-9
layer_with_weights-9
layer-10
layer_with_weights-10
layer-11
layer_with_weights-11
layer-12
layer_with_weights-12
layer-13
layer_with_weights-13
layer-14
layer_with_weights-14
layer-15
layer_with_weights-15
layer-16
layer_with_weights-16
layer-17
layer_with_weights-17
layer-18
layer_with_weights-18
layer-19
layer_with_weights-19
layer-20
	variables
trainable_variables
regularization_losses
	keras_api
__call__
*&call_and_return_all_conditional_losses
_default_save_signature
	optimizer
loss

signatures"
_tf_keras_network
"
_tf_keras_input_layer
�
 	keras_api
!
_keep_axis
"_reduce_axis
#_reduce_axis_mask
$_broadcast_shape
%mean
%
adapt_mean
&variance
&adapt_variance
	'count
(_adapt_function"
_tf_keras_layer
�
)	variables
*trainable_variables
+regularization_losses
,	keras_api
-__call__
*.&call_and_return_all_conditional_losses

/kernel
0bias"
_tf_keras_layer
�
1	variables
2trainable_variables
3regularization_losses
4	keras_api
5__call__
*6&call_and_return_all_conditional_losses

7kernel
8bias"
_tf_keras_layer
�
9	variables
:trainable_variables
;regularization_losses
<	keras_api
=__call__
*>&call_and_return_all_conditional_losses

?kernel
@bias"
_tf_keras_layer
�
A	variables
Btrainable_variables
Cregularization_losses
D	keras_api
E__call__
*F&call_and_return_all_conditional_losses

Gkernel
Hbias"
_tf_keras_layer
�
I	variables
Jtrainable_variables
Kregularization_losses
L	keras_api
M__call__
*N&call_and_return_all_conditional_losses

Okernel
Pbias"
_tf_keras_layer
�
Q	variables
Rtrainable_variables
Sregularization_losses
T	keras_api
U__call__
*V&call_and_return_all_conditional_losses

Wkernel
Xbias"
_tf_keras_layer
�
Y	variables
Ztrainable_variables
[regularization_losses
\	keras_api
]__call__
*^&call_and_return_all_conditional_losses

_kernel
`bias"
_tf_keras_layer
�
a	variables
btrainable_variables
cregularization_losses
d	keras_api
e__call__
*f&call_and_return_all_conditional_losses

gkernel
hbias"
_tf_keras_layer
�
i	variables
jtrainable_variables
kregularization_losses
l	keras_api
m__call__
*n&call_and_return_all_conditional_losses

okernel
pbias"
_tf_keras_layer
�
q	variables
rtrainable_variables
sregularization_losses
t	keras_api
u__call__
*v&call_and_return_all_conditional_losses

wkernel
xbias"
_tf_keras_layer
�
y	variables
ztrainable_variables
{regularization_losses
|	keras_api
}__call__
*~&call_and_return_all_conditional_losses

kernel
	�bias"
_tf_keras_layer
�
�	variables
�trainable_variables
�regularization_losses
�	keras_api
�__call__
+�&call_and_return_all_conditional_losses
�kernel
	�bias"
_tf_keras_layer
�
�	variables
�trainable_variables
�regularization_losses
�	keras_api
�__call__
+�&call_and_return_all_conditional_losses
�kernel
	�bias"
_tf_keras_layer
�
�	variables
�trainable_variables
�regularization_losses
�	keras_api
�__call__
+�&call_and_return_all_conditional_losses
�kernel
	�bias"
_tf_keras_layer
�
�	variables
�trainable_variables
�regularization_losses
�	keras_api
�__call__
+�&call_and_return_all_conditional_losses
�kernel
	�bias"
_tf_keras_layer
�
�	variables
�trainable_variables
�regularization_losses
�	keras_api
�__call__
+�&call_and_return_all_conditional_losses
�kernel
	�bias"
_tf_keras_layer
�
�	variables
�trainable_variables
�regularization_losses
�	keras_api
�__call__
+�&call_and_return_all_conditional_losses
�kernel
	�bias"
_tf_keras_layer
�
�	variables
�trainable_variables
�regularization_losses
�	keras_api
�__call__
+�&call_and_return_all_conditional_losses
�kernel
	�bias"
_tf_keras_layer
�
�	variables
�trainable_variables
�regularization_losses
�	keras_api
�__call__
+�&call_and_return_all_conditional_losses
�kernel
	�bias"
_tf_keras_layer
�
%0
&1
'2
/3
04
75
86
?7
@8
G9
H10
O11
P12
W13
X14
_15
`16
g17
h18
o19
p20
w21
x22
23
�24
�25
�26
�27
�28
�29
�30
�31
�32
�33
�34
�35
�36
�37
�38
�39
�40"
trackable_list_wrapper
�
/0
01
72
83
?4
@5
G6
H7
O8
P9
W10
X11
_12
`13
g14
h15
o16
p17
w18
x19
20
�21
�22
�23
�24
�25
�26
�27
�28
�29
�30
�31
�32
�33
�34
�35
�36
�37"
trackable_list_wrapper
 "
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
	variables
trainable_variables
regularization_losses
__call__
_default_save_signature
*&call_and_return_all_conditional_losses
&"call_and_return_conditional_losses"
_generic_user_object
�
�trace_0
�trace_1
�trace_2
�trace_32�
'__inference_model_layer_call_fn_2217294
'__inference_model_layer_call_fn_2218315
'__inference_model_layer_call_fn_2218406
'__inference_model_layer_call_fn_2217866�
���
FullArgSpec1
args)�&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0z�trace_1z�trace_2z�trace_3
�
�trace_0
�trace_1
�trace_2
�trace_32�
B__inference_model_layer_call_and_return_conditional_losses_2218553
B__inference_model_layer_call_and_return_conditional_losses_2218700
B__inference_model_layer_call_and_return_conditional_losses_2217975
B__inference_model_layer_call_and_return_conditional_losses_2218084�
���
FullArgSpec1
args)�&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0z�trace_1z�trace_2z�trace_3
�
�	capture_0
�	capture_1B�
"__inference__wrapped_model_2216864input"�
���
FullArgSpec
args� 
varargsjargs
varkwjkwargs
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�	capture_0z�	capture_1
�
�
_variables
�_iterations
�_learning_rate
�_index_dict
�
_momentums
�_velocities
�_update_step_xla"
experimentalOptimizer
 "
trackable_dict_wrapper
-
�serving_default"
signature_map
"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
:2mean
:2variance
:	 2count
�
�trace_02�
__inference_adapt_step_2218224�
���
FullArgSpec
args�

jiterator
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
.
/0
01"
trackable_list_wrapper
.
/0
01"
trackable_list_wrapper
 "
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
)	variables
*trainable_variables
+regularization_losses
-__call__
*.&call_and_return_all_conditional_losses
&."call_and_return_conditional_losses"
_generic_user_object
�
�trace_02�
$__inference_l2_layer_call_fn_2218709�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
�
�trace_02�
?__inference_l2_layer_call_and_return_conditional_losses_2218720�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
:@2	l2/kernel
:@2l2/bias
.
70
81"
trackable_list_wrapper
.
70
81"
trackable_list_wrapper
 "
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
1	variables
2trainable_variables
3regularization_losses
5__call__
*6&call_and_return_all_conditional_losses
&6"call_and_return_conditional_losses"
_generic_user_object
�
�trace_02�
$__inference_l3_layer_call_fn_2218729�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
�
�trace_02�
?__inference_l3_layer_call_and_return_conditional_losses_2218740�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
:@@2	l3/kernel
:@2l3/bias
.
?0
@1"
trackable_list_wrapper
.
?0
@1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
9	variables
:trainable_variables
;regularization_losses
=__call__
*>&call_and_return_all_conditional_losses
&>"call_and_return_conditional_losses"
_generic_user_object
�
�trace_02�
#__inference_x_layer_call_fn_2218749�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
�
�trace_02�
>__inference_x_layer_call_and_return_conditional_losses_2218760�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
:@@2x/kernel
:@2x/bias
.
G0
H1"
trackable_list_wrapper
.
G0
H1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
A	variables
Btrainable_variables
Cregularization_losses
E__call__
*F&call_and_return_all_conditional_losses
&F"call_and_return_conditional_losses"
_generic_user_object
�
�trace_02�
#__inference_1_layer_call_fn_2218769�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
�
�trace_02�
>__inference_1_layer_call_and_return_conditional_losses_2218780�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
:@@21/kernel
:@21/bias
.
O0
P1"
trackable_list_wrapper
.
O0
P1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
I	variables
Jtrainable_variables
Kregularization_losses
M__call__
*N&call_and_return_all_conditional_losses
&N"call_and_return_conditional_losses"
_generic_user_object
�
�trace_02�
#__inference_4_layer_call_fn_2218789�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
�
�trace_02�
>__inference_4_layer_call_and_return_conditional_losses_2218800�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
:@@24/kernel
:@24/bias
.
W0
X1"
trackable_list_wrapper
.
W0
X1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
Q	variables
Rtrainable_variables
Sregularization_losses
U__call__
*V&call_and_return_all_conditional_losses
&V"call_and_return_conditional_losses"
_generic_user_object
�
�trace_02�
#__inference_7_layer_call_fn_2218809�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
�
�trace_02�
>__inference_7_layer_call_and_return_conditional_losses_2218820�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
:@@27/kernel
:@27/bias
.
_0
`1"
trackable_list_wrapper
.
_0
`1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
Y	variables
Ztrainable_variables
[regularization_losses
]__call__
*^&call_and_return_all_conditional_losses
&^"call_and_return_conditional_losses"
_generic_user_object
�
�trace_02�
$__inference_10_layer_call_fn_2218829�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
�
�trace_02�
?__inference_10_layer_call_and_return_conditional_losses_2218840�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
:@@2	10/kernel
:@210/bias
.
g0
h1"
trackable_list_wrapper
.
g0
h1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
a	variables
btrainable_variables
cregularization_losses
e__call__
*f&call_and_return_all_conditional_losses
&f"call_and_return_conditional_losses"
_generic_user_object
�
�trace_02�
#__inference_2_layer_call_fn_2218849�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
�
�trace_02�
>__inference_2_layer_call_and_return_conditional_losses_2218860�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
:@@22/kernel
:@22/bias
.
o0
p1"
trackable_list_wrapper
.
o0
p1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
i	variables
jtrainable_variables
kregularization_losses
m__call__
*n&call_and_return_all_conditional_losses
&n"call_and_return_conditional_losses"
_generic_user_object
�
�trace_02�
#__inference_5_layer_call_fn_2218869�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
�
�trace_02�
>__inference_5_layer_call_and_return_conditional_losses_2218880�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
:@@25/kernel
:@25/bias
.
w0
x1"
trackable_list_wrapper
.
w0
x1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
q	variables
rtrainable_variables
sregularization_losses
u__call__
*v&call_and_return_all_conditional_losses
&v"call_and_return_conditional_losses"
_generic_user_object
�
�trace_02�
#__inference_8_layer_call_fn_2218889�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
�
�trace_02�
>__inference_8_layer_call_and_return_conditional_losses_2218900�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
:@@28/kernel
:@28/bias
/
0
�1"
trackable_list_wrapper
/
0
�1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
y	variables
ztrainable_variables
{regularization_losses
}__call__
*~&call_and_return_all_conditional_losses
&~"call_and_return_conditional_losses"
_generic_user_object
�
�trace_02�
$__inference_11_layer_call_fn_2218909�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
�
�trace_02�
?__inference_11_layer_call_and_return_conditional_losses_2218920�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
:@@2	11/kernel
:@211/bias
0
�0
�1"
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
�	variables
�trainable_variables
�regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
�
�trace_02�
#__inference_3_layer_call_fn_2218929�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
�
�trace_02�
>__inference_3_layer_call_and_return_conditional_losses_2218940�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
:@ 23/kernel
: 23/bias
0
�0
�1"
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
�	variables
�trainable_variables
�regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
�
�trace_02�
#__inference_6_layer_call_fn_2218949�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
�
�trace_02�
>__inference_6_layer_call_and_return_conditional_losses_2218960�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
:@ 26/kernel
: 26/bias
0
�0
�1"
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
�	variables
�trainable_variables
�regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
�
�trace_02�
#__inference_9_layer_call_fn_2218969�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
�
�trace_02�
>__inference_9_layer_call_and_return_conditional_losses_2218980�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
:@ 29/kernel
: 29/bias
0
�0
�1"
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
�	variables
�trainable_variables
�regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
�
�trace_02�
$__inference_12_layer_call_fn_2218989�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
�
�trace_02�
?__inference_12_layer_call_and_return_conditional_losses_2219000�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
:@ 2	12/kernel
: 212/bias
0
�0
�1"
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
�	variables
�trainable_variables
�regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
�
�trace_02�
#__inference_T_layer_call_fn_2219009�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
�
�trace_02�
>__inference_T_layer_call_and_return_conditional_losses_2219020�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
: 2T/kernel
:2T/bias
0
�0
�1"
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
�	variables
�trainable_variables
�regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
�
�trace_02�
#__inference_P_layer_call_fn_2219029�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
�
�trace_02�
>__inference_P_layer_call_and_return_conditional_losses_2219040�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
: 2P/kernel
:2P/bias
0
�0
�1"
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
�	variables
�trainable_variables
�regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
�
�trace_02�
%__inference_phi_layer_call_fn_2219049�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
�
�trace_02�
@__inference_phi_layer_call_and_return_conditional_losses_2219060�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
: 2
phi/kernel
:2phi/bias
0
�0
�1"
trackable_list_wrapper
0
�0
�1"
trackable_list_wrapper
 "
trackable_list_wrapper
�
�non_trainable_variables
�layers
�metrics
 �layer_regularization_losses
�layer_metrics
�	variables
�trainable_variables
�regularization_losses
�__call__
+�&call_and_return_all_conditional_losses
'�"call_and_return_conditional_losses"
_generic_user_object
�
�trace_02�
#__inference_c_layer_call_fn_2219069�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
�
�trace_02�
>__inference_c_layer_call_and_return_conditional_losses_2219080�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�trace_0
: 2c/kernel
:2c/bias
5
%0
&1
'2"
trackable_list_wrapper
�
0
1
2
3
4
5
6
7
	8

9
10
11
12
13
14
15
16
17
18
19
20"
trackable_list_wrapper
H
�0
�1
�2
�3
�4"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
�
�	capture_0
�	capture_1B�
'__inference_model_layer_call_fn_2217294input"�
���
FullArgSpec1
args)�&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�	capture_0z�	capture_1
�
�	capture_0
�	capture_1B�
'__inference_model_layer_call_fn_2218315inputs"�
���
FullArgSpec1
args)�&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�	capture_0z�	capture_1
�
�	capture_0
�	capture_1B�
'__inference_model_layer_call_fn_2218406inputs"�
���
FullArgSpec1
args)�&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�	capture_0z�	capture_1
�
�	capture_0
�	capture_1B�
'__inference_model_layer_call_fn_2217866input"�
���
FullArgSpec1
args)�&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�	capture_0z�	capture_1
�
�	capture_0
�	capture_1B�
B__inference_model_layer_call_and_return_conditional_losses_2218553inputs"�
���
FullArgSpec1
args)�&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�	capture_0z�	capture_1
�
�	capture_0
�	capture_1B�
B__inference_model_layer_call_and_return_conditional_losses_2218700inputs"�
���
FullArgSpec1
args)�&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�	capture_0z�	capture_1
�
�	capture_0
�	capture_1B�
B__inference_model_layer_call_and_return_conditional_losses_2217975input"�
���
FullArgSpec1
args)�&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�	capture_0z�	capture_1
�
�	capture_0
�	capture_1B�
B__inference_model_layer_call_and_return_conditional_losses_2218084input"�
���
FullArgSpec1
args)�&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults�
p 

 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�	capture_0z�	capture_1
!J	
Const_1jtf.TrackableConstant
J
Constjtf.TrackableConstant
�
�0
�1
�2
�3
�4
�5
�6
�7
�8
�9
�10
�11
�12
�13
�14
�15
�16
�17
�18
�19
�20
�21
�22
�23
�24
�25
�26
�27
�28
�29
�30
�31
�32
�33
�34
�35
�36
�37
�38
�39
�40
�41
�42
�43
�44
�45
�46
�47
�48
�49
�50
�51
�52
�53
�54
�55
�56
�57
�58
�59
�60
�61
�62
�63
�64
�65
�66
�67
�68
�69
�70
�71
�72
�73
�74
�75
�76"
trackable_list_wrapper
:	 2	iteration
: 2learning_rate
 "
trackable_dict_wrapper
�
�0
�1
�2
�3
�4
�5
�6
�7
�8
�9
�10
�11
�12
�13
�14
�15
�16
�17
�18
�19
�20
�21
�22
�23
�24
�25
�26
�27
�28
�29
�30
�31
�32
�33
�34
�35
�36
�37"
trackable_list_wrapper
�
�0
�1
�2
�3
�4
�5
�6
�7
�8
�9
�10
�11
�12
�13
�14
�15
�16
�17
�18
�19
�20
�21
�22
�23
�24
�25
�26
�27
�28
�29
�30
�31
�32
�33
�34
�35
�36
�37"
trackable_list_wrapper
�2��
���
FullArgSpec2
args*�'
jself

jgradient

jvariable
jkey
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 0
�
�	capture_0
�	capture_1B�
%__inference_signature_wrapper_2218179input"�
���
FullArgSpec
args� 
varargs
 
varkwjkwargs
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 z�	capture_0z�	capture_1
�B�
__inference_adapt_step_2218224iterator"�
���
FullArgSpec
args�

jiterator
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
$__inference_l2_layer_call_fn_2218709inputs"�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
?__inference_l2_layer_call_and_return_conditional_losses_2218720inputs"�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
$__inference_l3_layer_call_fn_2218729inputs"�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
?__inference_l3_layer_call_and_return_conditional_losses_2218740inputs"�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
#__inference_x_layer_call_fn_2218749inputs"�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
>__inference_x_layer_call_and_return_conditional_losses_2218760inputs"�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
#__inference_1_layer_call_fn_2218769inputs"�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
>__inference_1_layer_call_and_return_conditional_losses_2218780inputs"�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
#__inference_4_layer_call_fn_2218789inputs"�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
>__inference_4_layer_call_and_return_conditional_losses_2218800inputs"�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
#__inference_7_layer_call_fn_2218809inputs"�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
>__inference_7_layer_call_and_return_conditional_losses_2218820inputs"�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
$__inference_10_layer_call_fn_2218829inputs"�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
?__inference_10_layer_call_and_return_conditional_losses_2218840inputs"�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
#__inference_2_layer_call_fn_2218849inputs"�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
>__inference_2_layer_call_and_return_conditional_losses_2218860inputs"�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
#__inference_5_layer_call_fn_2218869inputs"�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
>__inference_5_layer_call_and_return_conditional_losses_2218880inputs"�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
#__inference_8_layer_call_fn_2218889inputs"�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
>__inference_8_layer_call_and_return_conditional_losses_2218900inputs"�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
$__inference_11_layer_call_fn_2218909inputs"�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
?__inference_11_layer_call_and_return_conditional_losses_2218920inputs"�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
#__inference_3_layer_call_fn_2218929inputs"�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
>__inference_3_layer_call_and_return_conditional_losses_2218940inputs"�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
#__inference_6_layer_call_fn_2218949inputs"�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
>__inference_6_layer_call_and_return_conditional_losses_2218960inputs"�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
#__inference_9_layer_call_fn_2218969inputs"�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
>__inference_9_layer_call_and_return_conditional_losses_2218980inputs"�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
$__inference_12_layer_call_fn_2218989inputs"�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
?__inference_12_layer_call_and_return_conditional_losses_2219000inputs"�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
#__inference_T_layer_call_fn_2219009inputs"�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
>__inference_T_layer_call_and_return_conditional_losses_2219020inputs"�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
#__inference_P_layer_call_fn_2219029inputs"�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
>__inference_P_layer_call_and_return_conditional_losses_2219040inputs"�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
%__inference_phi_layer_call_fn_2219049inputs"�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
@__inference_phi_layer_call_and_return_conditional_losses_2219060inputs"�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
�B�
#__inference_c_layer_call_fn_2219069inputs"�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
�B�
>__inference_c_layer_call_and_return_conditional_losses_2219080inputs"�
���
FullArgSpec
args�
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs� 
kwonlydefaults
 
annotations� *
 
R
�	variables
�	keras_api

�total

�count"
_tf_keras_metric
R
�	variables
�	keras_api

�total

�count"
_tf_keras_metric
R
�	variables
�	keras_api

�total

�count"
_tf_keras_metric
R
�	variables
�	keras_api

�total

�count"
_tf_keras_metric
R
�	variables
�	keras_api

�total

�count"
_tf_keras_metric
 :@2Adam/m/l2/kernel
 :@2Adam/v/l2/kernel
:@2Adam/m/l2/bias
:@2Adam/v/l2/bias
 :@@2Adam/m/l3/kernel
 :@@2Adam/v/l3/kernel
:@2Adam/m/l3/bias
:@2Adam/v/l3/bias
:@@2Adam/m/x/kernel
:@@2Adam/v/x/kernel
:@2Adam/m/x/bias
:@2Adam/v/x/bias
:@@2Adam/m/1/kernel
:@@2Adam/v/1/kernel
:@2Adam/m/1/bias
:@2Adam/v/1/bias
:@@2Adam/m/4/kernel
:@@2Adam/v/4/kernel
:@2Adam/m/4/bias
:@2Adam/v/4/bias
:@@2Adam/m/7/kernel
:@@2Adam/v/7/kernel
:@2Adam/m/7/bias
:@2Adam/v/7/bias
 :@@2Adam/m/10/kernel
 :@@2Adam/v/10/kernel
:@2Adam/m/10/bias
:@2Adam/v/10/bias
:@@2Adam/m/2/kernel
:@@2Adam/v/2/kernel
:@2Adam/m/2/bias
:@2Adam/v/2/bias
:@@2Adam/m/5/kernel
:@@2Adam/v/5/kernel
:@2Adam/m/5/bias
:@2Adam/v/5/bias
:@@2Adam/m/8/kernel
:@@2Adam/v/8/kernel
:@2Adam/m/8/bias
:@2Adam/v/8/bias
 :@@2Adam/m/11/kernel
 :@@2Adam/v/11/kernel
:@2Adam/m/11/bias
:@2Adam/v/11/bias
:@ 2Adam/m/3/kernel
:@ 2Adam/v/3/kernel
: 2Adam/m/3/bias
: 2Adam/v/3/bias
:@ 2Adam/m/6/kernel
:@ 2Adam/v/6/kernel
: 2Adam/m/6/bias
: 2Adam/v/6/bias
:@ 2Adam/m/9/kernel
:@ 2Adam/v/9/kernel
: 2Adam/m/9/bias
: 2Adam/v/9/bias
 :@ 2Adam/m/12/kernel
 :@ 2Adam/v/12/kernel
: 2Adam/m/12/bias
: 2Adam/v/12/bias
: 2Adam/m/T/kernel
: 2Adam/v/T/kernel
:2Adam/m/T/bias
:2Adam/v/T/bias
: 2Adam/m/P/kernel
: 2Adam/v/P/kernel
:2Adam/m/P/bias
:2Adam/v/P/bias
!: 2Adam/m/phi/kernel
!: 2Adam/v/phi/kernel
:2Adam/m/phi/bias
:2Adam/v/phi/bias
: 2Adam/m/c/kernel
: 2Adam/v/c/kernel
:2Adam/m/c/bias
:2Adam/v/c/bias
0
�0
�1"
trackable_list_wrapper
.
�	variables"
_generic_user_object
:  (2total
:  (2count
0
�0
�1"
trackable_list_wrapper
.
�	variables"
_generic_user_object
:  (2total
:  (2count
0
�0
�1"
trackable_list_wrapper
.
�	variables"
_generic_user_object
:  (2total
:  (2count
0
�0
�1"
trackable_list_wrapper
.
�	variables"
_generic_user_object
:  (2total
:  (2count
0
�0
�1"
trackable_list_wrapper
.
�	variables"
_generic_user_object
:  (2total
:  (2count�
?__inference_10_layer_call_and_return_conditional_losses_2218840c_`/�,
%�"
 �
inputs���������@
� ",�)
"�
tensor_0���������@
� �
$__inference_10_layer_call_fn_2218829X_`/�,
%�"
 �
inputs���������@
� "!�
unknown���������@�
?__inference_11_layer_call_and_return_conditional_losses_2218920d�/�,
%�"
 �
inputs���������@
� ",�)
"�
tensor_0���������@
� �
$__inference_11_layer_call_fn_2218909Y�/�,
%�"
 �
inputs���������@
� "!�
unknown���������@�
?__inference_12_layer_call_and_return_conditional_losses_2219000e��/�,
%�"
 �
inputs���������@
� ",�)
"�
tensor_0��������� 
� �
$__inference_12_layer_call_fn_2218989Z��/�,
%�"
 �
inputs���������@
� "!�
unknown��������� �
>__inference_1_layer_call_and_return_conditional_losses_2218780cGH/�,
%�"
 �
inputs���������@
� ",�)
"�
tensor_0���������@
� 
#__inference_1_layer_call_fn_2218769XGH/�,
%�"
 �
inputs���������@
� "!�
unknown���������@�
>__inference_2_layer_call_and_return_conditional_losses_2218860cgh/�,
%�"
 �
inputs���������@
� ",�)
"�
tensor_0���������@
� 
#__inference_2_layer_call_fn_2218849Xgh/�,
%�"
 �
inputs���������@
� "!�
unknown���������@�
>__inference_3_layer_call_and_return_conditional_losses_2218940e��/�,
%�"
 �
inputs���������@
� ",�)
"�
tensor_0��������� 
� �
#__inference_3_layer_call_fn_2218929Z��/�,
%�"
 �
inputs���������@
� "!�
unknown��������� �
>__inference_4_layer_call_and_return_conditional_losses_2218800cOP/�,
%�"
 �
inputs���������@
� ",�)
"�
tensor_0���������@
� 
#__inference_4_layer_call_fn_2218789XOP/�,
%�"
 �
inputs���������@
� "!�
unknown���������@�
>__inference_5_layer_call_and_return_conditional_losses_2218880cop/�,
%�"
 �
inputs���������@
� ",�)
"�
tensor_0���������@
� 
#__inference_5_layer_call_fn_2218869Xop/�,
%�"
 �
inputs���������@
� "!�
unknown���������@�
>__inference_6_layer_call_and_return_conditional_losses_2218960e��/�,
%�"
 �
inputs���������@
� ",�)
"�
tensor_0��������� 
� �
#__inference_6_layer_call_fn_2218949Z��/�,
%�"
 �
inputs���������@
� "!�
unknown��������� �
>__inference_7_layer_call_and_return_conditional_losses_2218820cWX/�,
%�"
 �
inputs���������@
� ",�)
"�
tensor_0���������@
� 
#__inference_7_layer_call_fn_2218809XWX/�,
%�"
 �
inputs���������@
� "!�
unknown���������@�
>__inference_8_layer_call_and_return_conditional_losses_2218900cwx/�,
%�"
 �
inputs���������@
� ",�)
"�
tensor_0���������@
� 
#__inference_8_layer_call_fn_2218889Xwx/�,
%�"
 �
inputs���������@
� "!�
unknown���������@�
>__inference_9_layer_call_and_return_conditional_losses_2218980e��/�,
%�"
 �
inputs���������@
� ",�)
"�
tensor_0��������� 
� �
#__inference_9_layer_call_fn_2218969Z��/�,
%�"
 �
inputs���������@
� "!�
unknown��������� �
>__inference_P_layer_call_and_return_conditional_losses_2219040e��/�,
%�"
 �
inputs��������� 
� ",�)
"�
tensor_0���������
� �
#__inference_P_layer_call_fn_2219029Z��/�,
%�"
 �
inputs��������� 
� "!�
unknown����������
>__inference_T_layer_call_and_return_conditional_losses_2219020e��/�,
%�"
 �
inputs��������� 
� ",�)
"�
tensor_0���������
� �
#__inference_T_layer_call_fn_2219009Z��/�,
%�"
 �
inputs��������� 
� "!�
unknown����������
"__inference__wrapped_model_2216864�;��/078?@_`WXOPGH�wxopgh����������������.�+
$�!
�
input���������
� "���
 
P�
p���������
 
T�
t���������
 
c�
c���������
$
phi�
phi���������p
__inference_adapt_step_2218224N'%&C�@
9�6
4�1�
����������IteratorSpec 
� "
 �
>__inference_c_layer_call_and_return_conditional_losses_2219080e��/�,
%�"
 �
inputs��������� 
� ",�)
"�
tensor_0���������
� �
#__inference_c_layer_call_fn_2219069Z��/�,
%�"
 �
inputs��������� 
� "!�
unknown����������
?__inference_l2_layer_call_and_return_conditional_losses_2218720c/0/�,
%�"
 �
inputs���������
� ",�)
"�
tensor_0���������@
� �
$__inference_l2_layer_call_fn_2218709X/0/�,
%�"
 �
inputs���������
� "!�
unknown���������@�
?__inference_l3_layer_call_and_return_conditional_losses_2218740c78/�,
%�"
 �
inputs���������@
� ",�)
"�
tensor_0���������@
� �
$__inference_l3_layer_call_fn_2218729X78/�,
%�"
 �
inputs���������@
� "!�
unknown���������@�
B__inference_model_layer_call_and_return_conditional_losses_2217975�;��/078?@_`WXOPGH�wxopgh����������������6�3
,�)
�
input���������
p 

 
� "���
���
$�!

tensor_0_0���������
$�!

tensor_0_1���������
$�!

tensor_0_2���������
$�!

tensor_0_3���������
� �
B__inference_model_layer_call_and_return_conditional_losses_2218084�;��/078?@_`WXOPGH�wxopgh����������������6�3
,�)
�
input���������
p

 
� "���
���
$�!

tensor_0_0���������
$�!

tensor_0_1���������
$�!

tensor_0_2���������
$�!

tensor_0_3���������
� �
B__inference_model_layer_call_and_return_conditional_losses_2218553�;��/078?@_`WXOPGH�wxopgh����������������7�4
-�*
 �
inputs���������
p 

 
� "���
���
$�!

tensor_0_0���������
$�!

tensor_0_1���������
$�!

tensor_0_2���������
$�!

tensor_0_3���������
� �
B__inference_model_layer_call_and_return_conditional_losses_2218700�;��/078?@_`WXOPGH�wxopgh����������������7�4
-�*
 �
inputs���������
p

 
� "���
���
$�!

tensor_0_0���������
$�!

tensor_0_1���������
$�!

tensor_0_2���������
$�!

tensor_0_3���������
� �
'__inference_model_layer_call_fn_2217294�;��/078?@_`WXOPGH�wxopgh����������������6�3
,�)
�
input���������
p 

 
� "���
"�
tensor_0���������
"�
tensor_1���������
"�
tensor_2���������
"�
tensor_3����������
'__inference_model_layer_call_fn_2217866�;��/078?@_`WXOPGH�wxopgh����������������6�3
,�)
�
input���������
p

 
� "���
"�
tensor_0���������
"�
tensor_1���������
"�
tensor_2���������
"�
tensor_3����������
'__inference_model_layer_call_fn_2218315�;��/078?@_`WXOPGH�wxopgh����������������7�4
-�*
 �
inputs���������
p 

 
� "���
"�
tensor_0���������
"�
tensor_1���������
"�
tensor_2���������
"�
tensor_3����������
'__inference_model_layer_call_fn_2218406�;��/078?@_`WXOPGH�wxopgh����������������7�4
-�*
 �
inputs���������
p

 
� "���
"�
tensor_0���������
"�
tensor_1���������
"�
tensor_2���������
"�
tensor_3����������
@__inference_phi_layer_call_and_return_conditional_losses_2219060e��/�,
%�"
 �
inputs��������� 
� ",�)
"�
tensor_0���������
� �
%__inference_phi_layer_call_fn_2219049Z��/�,
%�"
 �
inputs��������� 
� "!�
unknown����������
%__inference_signature_wrapper_2218179�;��/078?@_`WXOPGH�wxopgh����������������7�4
� 
-�*
(
input�
input���������"���
 
P�
p���������
 
T�
t���������
 
c�
c���������
$
phi�
phi����������
>__inference_x_layer_call_and_return_conditional_losses_2218760c?@/�,
%�"
 �
inputs���������@
� ",�)
"�
tensor_0���������@
� 
#__inference_x_layer_call_fn_2218749X?@/�,
%�"
 �
inputs���������@
� "!�
unknown���������@