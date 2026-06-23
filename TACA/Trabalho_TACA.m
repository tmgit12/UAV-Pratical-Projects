%% Simulates OFDM
EN=[0:2:22]'+0*100; en = 10 .^(EN/10) ;
N=512;
NSlot=3255; % para fazer aproximadamene 10e6 bits com 3072 bits/slot precisamos de 10e6/3072=3255 slots
CHANNEL='AWGN';
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
    % --- Huffman Stream to 64-QAM Mapping ---
    bit_idx_start = (nn-1)*bits_per_slot + 1;
    bit_idx_end = nn*bits_per_slot;
    bits_Tx = full_encoded_stream(bit_idx_start:bit_idx_end).';
    
    % Modulate in 64-QAM
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


%% Parte 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configuração
CODING_SCHEME = 'CONV';      % 'UNC' (Uncoded), 'HAMMING', 'CONV', 'CONCAT', 'TURBO'
CHANNEL       = 'RAYL';      % 'AWGN' ou 'RAYL'
ESTIMATION    = 'IMPERFECT'; % 'PERFECT' ou 'IMPERFECT'
M             = 64;          % 4 para QPSK, 64 para 64-QAM
EN_dB         = 0:2:30;      % Eb/No em dB (até 30dB para Rayleigh/Rice, até 7 para AWGN)
NSlot         = 500;

% Parametros OFDM
N=512;       % Subcarriers
Ts=4e-6;      % Block duration
Tg=0.2 * Ts;  % Cyclic prefix durration
L=1; % L-th order diversity
k_mod = log2(M);% Bits per symbol


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Configurar FEC e baralhador
info_frame_len = 1122; % Base information bits per frame (multiple of 11 for Hamming and of 6 for for 64 QAM)

switch CODING_SCHEME
    case 'UNC'
        code_rate = 1;
        enc_len = info_frame_len;
        
    case 'HAMMING'
        code_rate = 11/15;
        % O matlab já trazia um codificador de Hamming
        hHamEnc = comm.HammingEncoder(15, 11);
        hHamDec = comm.HammingDecoder(15, 11);
        enc_len = info_frame_len * (15/11);
        
    case 'CONV'
        code_rate = 1/2;
        trellis = poly2trellis(7, [171 133]);
        enc_len = info_frame_len * 2;
        
    case 'CONCAT'
        % Junção do código de blocos com o convolucional
        code_rate = (11/15) * (1/2);
        hHamEnc = comm.HammingEncoder(15, 11);
        hHamDec = comm.HammingDecoder(15, 11);
        trellis = poly2trellis(7, [171 133]);
        enc_len = info_frame_len * (15/11) * 2;
        
    case 'TURBO'
        code_rate = 1/3;
        turboTrellis = poly2trellis(4, [13 15 17], 13);
        intrlvrInd = randperm(RandStream('mt19937ar','Seed',11), info_frame_len);
        hTEnc = comm.TurboEncoder('TrellisStructure', turboTrellis, 'InterleaverIndices', intrlvrInd);
        hTDec = comm.TurboDecoder('TrellisStructure', turboTrellis, 'InterleaverIndices', intrlvrInd, 'NumIterations', 4);
        enc_len = info_frame_len * 3 + 12; % Turbo adds 12 tail bits
end

% Baralhador
rng(12345); % Fixed seed for interleaver
permVec = randperm(enc_len)';


% Numero de bits necessarios, repeticao da mensagem, os codifiadores
% recebem colunas
total_info_bits_needed = NSlot * info_frame_len;
num_repeats = ceil(total_info_bits_needed / length(base_encoded_msg));
phase2_input_stream = repmat(base_encoded_msg, 1, num_repeats);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulação
en = 10 .^ (EN_dB / 10);
NEN = length(EN_dB);
BER = zeros(NEN, 1);
FER = zeros(NEN, 1); % Frame Error Rate = Retransmission Probability

fprintf('Starting %s simulation over %s channel (%s estimation)\n', CODING_SCHEME, CHANNEL, ESTIMATION);

for nEN = 1:NEN
    
    % Eb=1 QPSK; Eb=2.5 16QAM; Eb=7 64QAM
    % É NECESSÁRIO AJUSTAR E_b COM BASE NA MODULAÇÃO
    Eb = 1 / (k_mod * code_rate);
    sigma = sqrt(Eb / 2 / en(nEN));
    current_noise_var = 2 * (sigma^2);
    
    bit_err_count = 0;
    frame_err_count = 0;
    total_bits = 0;
    
    for nn = 1:NSlot
        
        % Gerar bits e codificar
        start_idx = (nn-1) * info_frame_len + 1;
        end_idx = nn * info_frame_len;
        info_bits = phase2_input_stream(start_idx:end_idx);
        info_bits = info_bits(:);
        
        switch CODING_SCHEME
            case 'UNC'
                enc_bits = info_bits;
            case 'HAMMING'
                enc_bits = step(hHamEnc, info_bits);
            case 'CONV'
                enc_bits = convenc(info_bits, trellis);
            case 'CONCAT'
                outer_bits = step(hHamEnc, info_bits);
                enc_bits = convenc(outer_bits, trellis);
            case 'TURBO'
                enc_bits = step(hTEnc, info_bits);
        end
        
        % baralhar e modular
        int_bits = intrlv(enc_bits(:), permVec);
        tx_sym = qammod(int_bits(:), M, 'InputType', 'bit', 'UnitAveragePower', true);
        
        % juntar em blocos OFDM
        num_syms = length(tx_sym);
        num_ofdm_blocks = ceil(num_syms / N);
        pad_len = num_ofdm_blocks * N - num_syms;
        tx_sym_padded = [tx_sym(:); zeros(pad_len, 1)]; % Pad to fit OFDM subcarriers
        
        rx_sym_padded = zeros(size(tx_sym_padded));
        
        % transmissão OFDM
        for b = 1:num_ofdm_blocks
            idx = (b-1)*N + 1 : b*N;
            Ak_Tx = tx_sym_padded(idx);
            
            % Channel Generation
            if strcmp(CHANNEL, 'RAYL')
                Hk = (randn(N,L) + 1i*randn(N,L)) / sqrt(2);
            else % AWGN
                Hk = ones(N,L) .* exp(1i*2*pi*rand(N,L));
            end
            sH2k = sum(abs(Hk).^2, 2);
            
            % Channel Estimation
            if strcmp(ESTIMATION, 'IMPERFECT')
                amp_err = -0.1 + 0.2 * rand(N, L);         % Max 10% amp error
                phase_err = deg2rad(-5 + 10 * rand(N, L)); % Max 5 deg phase error
                H_est = Hk .* (1 + amp_err) .* exp(1i * phase_err);
            else
                H_est = Hk;
            end
            sH2k_est = sum(abs(H_est).^2, 2);
            
            % Additive Noise & Equalization
            Yk = zeros(N, L);
            YIk = zeros(N, 1);
            for l = 1:L
                Yk(:,l) = Ak_Tx .* Hk(:,l) + (randn(N,1) + 1i*randn(N,1)) * sigma;
                YIk = YIk + Yk(:,l) .* conj(H_est(:,l)); % Equalize with estimate
            end
            YIk = YIk ./ sH2k_est;
            rx_sym_padded(idx) = YIk;
        end
        
        % desmodular
        rx_sym = rx_sym_padded(1:num_syms); % Remove padding
        
        if strcmp(CODING_SCHEME, 'TURBO')
            % Turbo needs Soft LLRs
            demod_out = qamdemod(rx_sym(:), M, 'OutputType', 'approxllr', ...
                'UnitAveragePower', true, 'NoiseVariance', current_noise_var);
        else
            % Others use Hard Bits for simplicity
            demod_out = qamdemod(rx_sym(:), M, 'OutputType', 'bit', ...
                'UnitAveragePower', true);
        end
        
        % desbaralhar e descodificar
        deint_out = deintrlv(demod_out(:), permVec);
        
        switch CODING_SCHEME
            case 'UNC'
                dec_bits = deint_out;
            case 'HAMMING'
                dec_bits = step(hHamDec, deint_out);
            case 'CONV'
                dec_bits = vitdec(deint_out, trellis, 35, 'trunc', 'hard');
            case 'CONCAT'
                inner_dec = vitdec(deint_out, trellis, 35, 'trunc', 'hard');
                dec_bits = step(hHamDec, inner_dec);
            case 'TURBO'
                dec_bits = step(hTDec, deint_out);
        end
        
        % erros
        errs_in_frame = sum(info_bits(:) ~= dec_bits(:));
        bit_err_count = bit_err_count + errs_in_frame;
        total_bits = total_bits + info_frame_len;
        
        if errs_in_frame > 0
            frame_err_count = frame_err_count + 1;
        end
    end
    
    BER(nEN) = bit_err_count / total_bits;
    FER(nEN) = frame_err_count / NSlot; % Retransmission probability
    
    fprintf('Eb/No: %2d dB | BER: %e | FER (Retransmissions): %e\n', EN_dB(nEN), BER(nEN), FER(nEN));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resultados visuais
figure;
subplot(1,2,1);
semilogy(EN_dB, BER, 'LineWidth', 2, 'Marker', 'o');
grid on; xlabel('E_b/N_0 (dB)'); ylabel('BER');
title(sprintf('BER: %s | %s | %s Est', CODING_SCHEME, CHANNEL, ESTIMATION));
ylim([1e-5 1]);

subplot(1,2,2);
semilogy(EN_dB, FER, 'LineWidth', 2, 'Marker', 's', 'Color', 'r');
grid on; xlabel('E_b/N_0 (dB)'); ylabel('FER (Retransmission Prob)');
title('Retransmission Probability');
ylim([1e-5 1]);
