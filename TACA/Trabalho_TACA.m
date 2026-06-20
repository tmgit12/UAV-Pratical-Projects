%% Simulates OFDM
EN=[0:2:22]'+0*100; en = 10 .^(EN/10) ;
N=512;
NSlot=3255; % para fazer aproximadamene 10e6 bits com 3072 bits/slot precisamos de 10e6/3072=3255 slots
CHANNEL='RAYL';
L=1; % L-th order diversity
Ts=4e-6; % Block duration
Tg=0.2*Ts; % Cyclic prefix durration
f=[-N/2:N/2-1]'/Ts; % frequencies

%% Mensagem e Huffman
msg = ['convolucional codes are still used in telemetry channels in ' ...
    'aeronautical applications and deep space communications due to ' ...
    'good performance and simplicity'];

alphabet = unique(msg);

prob = zeros(1, length(alphabet));

for i = 1:length(alphabet)
    prob(i) = sum(msg == alphabet(i)) / length(msg);
end
% sum(probabilities) Confirma-se q dá 1

%huffman dict
%bits=huffmancod

dict = huffmandict(num2cell(alphabet), prob);
% É preciso fazer num2cell(alphabet) para passar de 'abc...'  para 
% {'a', 'b', 'c', ...}
base_encoded_msg = huffmanenco(num2cell(msg), dict);


% Ritmos (bits/simbolo)
% Huffman
% nº de simbolos no dict
num_symbols = size(dict, 1);

% comprimento de cada simbolo
huffman_lengths = zeros(1, num_symbols);
for i = 1:num_symbols
    
    % Ver a palavra codificada na 2a coluna
    codeword = dict{i, 2};
    huffman_lengths(i) = length(codeword);
end
huffman_rate = sum(prob .* huffman_lengths');  % Expected value

% Codigo equiprovavel
% Log2(27) o alfabeto tem 27 simbolos únicos (a contar com o espaço)
equiprobable_rate = ceil(log2(length(alphabet))); 

% ASCII usa 8 bits/mensagem
ascii_rate = 8; 

% Eficiencia de compressao
compression_gain = (1 - (huffman_rate / equiprobable_rate)) * 100;

% Bits para a simulação
M = 64;                 % 64-QAM
k = log2(M);            % 6 bits per symbol
bits_per_slot = N * k;  % 512 * 6 = 3072 bits per OFDM block
total_bits_needed = bits_per_slot * NSlot;

% Repeat encoded message until we have enough bits for all slots
num_repeats = ceil(total_bits_needed / length(base_encoded_msg));
full_encoded_stream = repmat(base_encoded_msg, 1, num_repeats);
% ---------------------------------------------------

if (CHANNEL=='AWGN')
    alfa_med=1;tau=0;NRay=1;
elseif (CHANNEL=='REAL')
%    load pdp.940
%    load PDP_ChA.dat,pdp=PDP_ChA;
    load pdp_hipc.dat,pdp=pdp_hipc;
    tau=pdp(:,3);
    alpha_med=10 .^((pdp(:,1))/20);
    alpha_med=alpha_med/sqrt(sum(alpha_med.^2));
    NRay=length(alpha_med);
%    NRay=10;tau=[0:1:NRay-1]'*0.1/NRay*Ts; % Tg=0.1Ts
elseif (CHANNEL=='XTAP')
    NRay=32;
    tau=[0:NRay-1]'*Ts/N;
    alpha_med=ones(NRay,1);
    alpha_med=alpha_med/sqrt(sum(alpha_med.^2));
end;

Eb=7; % Eb=1 QPSK; Eb=2.5 16QAM; Eb=7 64QAM
sigma=sqrt(Eb/2 ./en); 
NSR=1/2 ./(en); 
NEN=length(EN);
NErr=zeros(NEN,1);
for nn=1:NSlot
    
    %rand('state',nn*1234567); randn('state',nn*1234567);
    % This means the same channel for each slot

    if (CHANNEL=='REAL')
        Hk=zeros(N,L); 
        for l=1:L
            alpha=alpha_med.*(randn(NRay,1)+j*randn(NRay,1))/sqrt(2);
            for nRay=1:NRay
                Hk(:,l)=Hk(:,l)+alpha(nRay)*exp(-j*2*pi*f*tau(nRay));
            end;
        end; 
    elseif (CHANNEL=='XTAP')
        Hk=zeros(N,L); 
        for l=1:L
            alpha=alpha_med.*(randn(NRay,1)+j*randn(NRay,1))/sqrt(2);
            for nRay=1:NRay
                Hk(:,l)=Hk(:,l)+alpha(nRay)*exp(-j*2*pi*f*tau(nRay));
            end;
        end;
    elseif (CHANNEL=='RRND')
        Hk=zeros(N,L); 
        tau=rand(NRay,1)*Tg;
            for l=1:L
                alpha=ones(NRay,1).*(randn(NRay,1)+j*randn(NRay,1))/sqrt(2*NRay);
                for nRay=1:NRay
                    Hk(:,l)=Hk(:,l)+alpha(nRay)*exp(-j*2*pi*f*tau(nRay));
                end;
            end; 
    elseif (CHANNEL=='RAYL')
        Hk=(randn(N,L)+j*randn(N,L))/sqrt(2);
    elseif (CHANNEL=='AWGN')
        Hk=ones(N,L).*exp(j*2*pi*rand(N,L));
    end;
    H2k=abs(Hk).^2;
    if (L==1) sH2k=H2k; else sH2k=sum(H2k')'; end;
    % --- CHANGED: 3. Huffman Stream to 64-QAM Mapping ---
    % Grab the specific chunk of bits for this slot
    bit_idx_start = (nn-1)*bits_per_slot + 1;
    bit_idx_end = nn*bits_per_slot;
    bits_Tx = full_encoded_stream(bit_idx_start:bit_idx_end).';
    
    % Modulate using 64-QAM
    %Ak_Tx=qammod(bits,64,InputType='bit')
    Ak_Tx = qammod(bits_Tx, M, 'InputType', 'bit');
    %Ak_Tx=sign(randn(N,1))+j*sign(randn(N,1)); % +/-1 +/-j
    an_Tx=fftshift(ifft(fftshift(Ak_Tx)));

    for nEN=1:NEN
        Yk=zeros(N,L);
        for l=1:L
            Yk(:,l)=Ak_Tx.*Hk(:,l)+(randn(N,1)+j*randn(N,1))*sigma(nEN);
        end;
        YIk=0;
        for l=1:L
            YIk = YIk +Yk(:,l).*conj(Hk(:,l));
        end;
        YIk=YIk./sH2k;

        
        %64-QAM Demodulation & Bit Error Calculation
        %bits_rx=qamdemod(YIK,M,OutputType='bit')
        %aux=abs(bits-bits_rx)
        bits_Rx = qamdemod(YIk, M, 'OutputType', 'bit');
        
        % Direct bit-by-bit comparison
        

        %Ak_Rx=sign(real(YIk))+j*sign(imag(YIk));
        %aux = sum( abs(real(Ak_Tx)-real(Ak_Rx)) + ...
                %abs(imag(Ak_Tx)-imag(Ak_Rx)) ) / 2 ;
        aux = sum(bits_Tx ~= bits_Rx);
        NErr(nEN,1)=NErr(nEN,1)+aux;
    end;

    if (rem(nn,100)==0) nn, end;
end;

% BER in Rayleigh channel and L-branch diversity [Proakis]
aux=sqrt(en./(1+en));Pb_tr=0;
% for l=0:L-1
%     Pb_tr=Pb_tr+combin(L-1+l,l)*((1+aux)/2).^l;
% end;
% Pb_tr=Pb_tr.*channel((1-aux)/2).^L;

% BER
Pb = NErr / (NSlot * bits_per_slot);
%%
% Theoretical BER for 64-QAM in AWGN (approximation)
% Note: Using standard qfunc since q_x is a custom function
% BER in AWGN 
PbAWGN=q_x(sqrt(2*L*en));

%figure;
semilogy(EN,Pb,'g-*',EN,PbAWGN,'b:')
xlabel('E_b/N_0(dB)'),ylabel('BER')
axis([0 20 1e-4 1])
%pause,clf;
