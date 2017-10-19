%%%%%%%%%%%%%%%%%%%%%%%%%%  
%  
% Revision Code: 8zNrNv  
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%  


## OVERVIEW

Welcome to the Simon Fraser University Cognitive Science Laboratory modeling repository.   

This readme will explain the coding standards, and internal revision labels we sometimes use to quickly describe the state of a function/script.

Models are added as submodules in this repository, while most other scripts can be committed directly.


## STANDARDS

* Publication critical code must be reviewed and verified by people other than the original authors.
* Commits should try to follow best practices:
	* https://gist.github.com/adeekshith/cd4c95a064977cdc6c50
* Documentation should be written in Markdown (pretty), reStructuredText (pretty + programmatic + complicated), or Latex (unusual for basic documentation but obviously the most powerful)
	* https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet
	* https://github.com/jgm/pandoc - convert between formats  
	* https://thomas-cokelaer.info/tutorials/sphinx/rest_syntax.html - learn reStructuredText  
	
* Code should be presentable in publication.

	To make Matlab code look pretty in Latex:

	https://github.com/Jubobs/matlab-prettifier

	To try and auto-convert Matlab expressions to Latex:

	http://people.csail.mit.edu/pgbovine/mathviz/
	https://github.com/pgbovine/mathviz

	To convert Latex to Matlab, try:

	https://github.com/ProjectTiresias/Tiresias/

	To auto-convert Matlab figures to Latex:

	https://github.com/matlab2tikz/

	Model configurations should be set by external json files, e.g.

	https://github.com/fangq/jsonlab  
	
* See the feature request grammar below as well.  
	
## Repository structure/organization  

Goal: Move everything from lab SVN to Git.  

### Org chart

1. Public (Public. Reviewed/Verified)  
	* LAG-1 (subtree)  
	* Data (git-lfs)    
2. Analyses (Private)  
	* Gnarly (subtree)  
	* SQL tools  
	* Experiments  
	* Models  
		* LAG-1  
3. InProgress (Private)  
	* Analyses  
	* Modeling  
		* RLAttn
	* Experiments  
		* CovertLearning  
4. Modeling (Public)  
	* Common  
		* ThirdParty  
			* Cosivina (Submodule. Mercurial. git-hg) 
		* CSLab  
5. Gnarly  
6. LAG-1  
7. Skillcraft  
8. Data Visualization  
9. Wiki  


## GRAMMAR FOR FEATURE REQUESTS

Every function or code snippet committed to this repository contains a header like in this example:

%%%%%%%%%%%%%%%%%%%%%%%%%%  
%  
% Revision Code: 3a5c8w0r0v  
%	- needs a bit of input decimation  
%	- code is in need of function compression  
%	- wildcard notes  
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%  

RevisionParse :=  

	[
	 InputDecimation = "regexp:^([1-9]|N)[a]?(.*?)[\s]*$"  
	 MagicNumberDecimation = "regexp:^([1-9]|N)[b]?(.*?)[\s]*$"  
	 CodeRedundancy = "regexp:^([1-9]|N)[c]?(.*?)[\s]*$"  
	 Reviewed = "regexp:^([1-9]|N)[r]?(.*?)[\s]*$"  
	 Verified = "regexp:^([1-9]|N)[v]?(.*?)[\s]*$"  
	 Rereview = "regexp:^([1-9]|N)[R]?(.*?)[\s]*$"  
	 Reverify = "regexp:^([1-9]|N)[V]?(.*?)[\s]*$"  
	 Wildcard = "regexp:^([1-9]|N)[w]?(.*?)[\s]*$"  
	]  

	UpdatePriority = InputDecimation + MagicNumberDecimation ...  
		CodeRedundancy + Reviewed + Verified  


%% InputDecimation

Automated function that increases priority linearly with the number of inputs.

%% MagicNumberDecimation

Automated analysis that increases priority linearly with the number non variable digit declarations.

%% CodeRedundancy

Automated function that increases priority linearly with the amount of disparate but redundant code segments 

%% Reviewed

0 if unreviewed, 1 if reviewed, N if Not applicable

%% Verified

0 if unverified, 1 if verified, N if Not applicable

%% ReReview needed

1 if ReReview needed, 0 otherwise

%% ReVerify needed

1 if ReReview needed, 0 otherwise

%% Wildcard note

0-9 level of urgency specified in a header note
