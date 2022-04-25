function [H, F] = dtft(h,Fs)
%DTFT   calculate DTFT at N equally spaced frequencies (numbers of sampling to be 2^n)
%----
%   Usage:   [H, F] = dtft(h, Fs)
%
%      h : finite-length input vector, whose length is L
%      N : number of frequencies for evaluation over [-pi,pi)
%              ==> constraint: N >= L 
%      H : DTFT values (complex)
%      F : (2nd output) vector of freqs where DTFT is computed
%           Single-Sided Ban Spectrum

%---------------------------------------------------------------
% copyright 1994, by C.S. Burrus, J.H. McClellan, A.V. Oppenheim,
% T.W. Parks, R.W. Schafer, & H.W. Schussler.  For use with the book
% "Computer-Based Exercises for Signal Processing Using MATLAB"
% (Prentice-Hall, 1994).
%---------------------------------------------------------------

%N = fix(N);
L = 2*length(h);  h = h(:);  %<-- for vectors ONLY !!!
%if( N < L )
%   error('DTFT: # data samples cannot exceed # freq samples')
%end
N=2^nextpow2(L); % FFT needs numbers of sampling to be 2^n
H = fft( h, N );
%Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.
P2 = abs(H/N);
P1 = P2(1:N/2+1);
P1(2:end-1) = 2*P1(2:end-1);
F =Fs*(0:(N/2))/N;
H = P1;  

