clc
clear all
close all



%[x,fs]=audioread('E:\Final Project-competetion\project results\bugging.wav');
%[x,fs]=audioread('E:\Final Project-competetion\project results\got.wav');
[x,fs]=audioread('E:\Final Project-competetion\project results\titon_2_07_SNR5.wav');
%[x,fs]=audioread('E:\Final Project-competetion\project results\yifen_2_01_SNR5.wav');


[m,n]=size(x);

wlen=1024
hop=1024/4
nfft=1024

% length of the signal
xlen = length(x);

% form a periodic hamming window

win = hamming(wlen);




% stft matrix estimation and preallocation
rown = ceil((1+nfft)/2);            % calculate the total number of rows
coln = 1+fix((xlen-wlen)/hop);      % calculate the total number of columns
stft = zeros(rown, coln);           % form the stft matrix

% initialize the signal time segment index
indx = 0;

% perform STFT
for col = 1:coln
    % windowing
    xw = x(indx+1:indx+wlen).*(win);
    
    % FFT
    X = fft(xw, nfft);
    
    % update the stft matrix
    stft(:, col) = X(1:rown);
    
    % update the index
    indx = indx + hop;
end
stft1=abs(stft);
D=stft1;

  t = (wlen/2:hop:wlen/2+(coln-1)*hop)/fs;
f = (0:rown-1)*fs/nfft;

  S1 = 20*log10(D);

% plot the spectrogram
figure(1);
surf(t, f, S1);
shading interp;
axis tight;
box on;
view(0, 90);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
xlabel('Time, s');
ylabel('Frequency, Hz');
title('Amplitude spectrogram of the signal');


X=D;
[M N]=size(D);
[m n] = size(D);
[M, N] = size(X);
   


unobserved = isnan(X);
    X(unobserved) = 0;
   normX = norm(X,'fro');
   

    % default arguments
  
        lambda = 1 / sqrt(max(M,N));
  
        mu =0.002;
   
        tol = 1e-6;
    
        max_iter =500;
    
    
    % initial solution
    L = zeros(M, N);
    S1 = zeros(M, N);
    Y = zeros(M, N);
    for iter = (1:max_iter)                                                            
       
        L = Do(1/mu, X - S1 + (1/mu)*Y);
        S1 = So(lambda/mu, X - L + (1/mu)*Y);
       
        Z = X - L - S1;
        Z(unobserved) = 0; % skip missing values
        Y = Y + mu*Z;
        
        err = norm(Z, 'fro') / normX;
        if (iter == 1) || (mod(iter, 10) == 0) || (err < tol)
            fprintf(1, 'iter: %04d\terr: %f\trank(L): %d\tcard(S): %d\n', ...
                    iter, err, rank(L), nnz(S1(~unobserved)));
        end
        if (err < tol) break; 
        end
    end



        
  



S=abs(L);
   S2 = 20*log10(S);

% plot the spectrogram
figure(2);
surf(t, f, S2);
shading interp;
axis tight;
box on;
view(0, 90);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
xlabel('Time, s');
ylabel('Frequency, Hz');
title('Amplitude spectrogram of the signal');

 S=abs(S1);
   S2 = 20*log10(S);

% plot the spectrogram
figure(3);
surf(t, f, S2);
shading interp;
axis tight;
box on;
view(0, 90);
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
xlabel('Time, s');
ylabel('Frequency, Hz');


for i=1:m
    for j=1:n
        if(abs(L(i,j))>(abs(S1(i,j))))
            D3(i,j)=1;
        else
            D3(i,j)=0;
        end
    end
end
   E_hatq=stft.*D3;;
   g=ones(m,n);
     L=stft.*(g-D3);
     
%% Istft and recovery for 4 songs
     
%    [swapnzsdfhsdfh,t] = istft(E_hatq, wlen, hop, nfft, fs);
%  audiowrite('E:\Final Project-competetion\project results\output\buggin1inst.wav',swapnzsdfhsdfh,fs);
%     [swapnilphsagsw,t] = istft(L, wlen, hop, nfft, fs);
%  audiowrite('E:\Final Project-competetion\project results\output\buggin1voc.wav',swapnilphsagsw,fs);


%  [swapnzsdfhsdfh,t] = istft(E_hatq, wlen, hop, nfft, fs);
%  audiowrite('E:\Final Project-competetion\project results\output\got1inst.wav',swapnzsdfhsdfh,fs);
%     [swapnilphsagsw,t] = istft(L, wlen, hop, nfft, fs);
%  audiowrite('E:\Final Project-competetion\project results\output\got1voc.wav',swapnilphsagsw,fs);

 
  [swapnzsdfhsdfh,t] = istft(E_hatq, wlen, hop, nfft, fs);
 audiowrite('E:\Final Project-competetion\project results\output\titan1inst.wav',swapnzsdfhsdfh,fs);
    [swapnilphsagsw,t] = istft(L, wlen, hop, nfft, fs);
 audiowrite('E:\Final Project-competetion\project results\output\titan1voc.wav',swapnilphsagsw,fs);

 
%   [swapnzsdfhsdfh,t] = istft(E_hatq, wlen, hop, nfft, fs);
%  audiowrite('E:\Final Project-competetion\project results\output\yifen1inst.wav',swapnzsdfhsdfh,fs);
%     [swapnilphsagsw,t] = istft(L, wlen, hop, nfft, fs);
%  audiowrite('E:\Final Project-competetion\project results\output\yifen1voc.wav',swapnilphsagsw,fs);
