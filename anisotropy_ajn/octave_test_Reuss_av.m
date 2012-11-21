%CReuss = octave_test_Reuss_av()
%
%function [CReuss] = octave_test_Reuss_av()
% Tests the module function CIJ_Reuss_av using the matrix manipulations within
% Octave.

	N=5 ;
	
	S = zeros(6,6) ;
	
	for n=1:N
		C{n} = zeros(6,6) ;
		for i=1:6
			for j=i:6
				C{n}(i,j) = n*i ;
				C{n}(j,i) = C{n}(i,j) ;
			end
		end
		
		S = S + inv(C{n})./N ;
	end
	
	CReuss = inv(S)

%end