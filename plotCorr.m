  
function o = plotCorr(F)
        % clc;
        % clear all;
        % close all;

        % Sequence and other parameters
        % X: a set of M sequences of length N
        M         = 5;
        N         = 100;

        % other parameters
        r         = 25;          % R-th root of unity

        normed    = 0;           % 0: un-normalized axis, 
                                 % 1: normalized axis
        ch_X      = 1;           % 1: random complex sequence a + ib (a,b are in [-1,1])
                                 % 2: binary sequence [1,-1]
                                 % 3: all ones 
                                 % 4: finite alphabate: rth roots of unity
                                 % 5: single r-th root ZadoffChu sequence
                                 % 6: chirp signal

        % Single sequence properties                         
        params    = [];
        params.M  = M;
        params.N  = N;
        params.R  = r;
        idx1      = 1;           % which pair you want to see the correlations
        idx2      = 1;
        if (idx1 > M || idx2 > M)
            error('sequence indices cannot be larger than M');
        end

        X = genSignal(ch_X, params);


        % cross-correlations
        [R, kr]  = correlation(F, 'a');
        % [C, kc]  = correlation(X, 'p');

        % plots
        if normed == 1
            temp_d   = (norm(X(:,idx1))*norm(X(:,idx2)));
            temp_str = 'Normalized amplitude'; 
        elseif normed == 0
            temp_d   = 1;
            temp_str = 'Un-normalized amplitude'; 
        end
        idx1 = 1
        b=1 ;
        dic = [ 'A' ' ' 'B' ' ' 'C' ' ' 'D' ] ;
        for idx1 = 1:2:7
            a=1 ;
            sum = 0 ;
            figure(); 
                for idx2 = 1:2:7

                    subplot(2,4,a); 
                                    stem(kr,squeeze( R(idx1,idx2,:) ) + squeeze( R(idx1+1,idx2+1,:) )    );
                                    title(sprintf('Aperiodic cross-correlation %s and %s', dic(idx1) , dic(idx2)) );
                                    xlabel('k'); ylabel(temp_str); grid on;
                    subplot(2,4,a+1); 
                                    stem(kr, squeeze( R(idx2,idx1,:) ) + squeeze( R(idx2+1,idx1+1,:) )  );
                                    title(sprintf('Aperiodic cross-correlation %s and %s ', dic(idx2), dic(idx1) ));
                                    xlabel('k'); ylabel(temp_str); grid on;
                    a = a + 2 ;
                    
%                     if idx1 ~= idx2 
%                         sum = sum + R(idx1,idx2,:)  +  R(idx1+1,idx2+1,:) ;
%                     end
                end 
%                 subplot(3,4,6+b);
%                 stem( kr , squeeze( sum ) ) ;

                suptitle( dic(idx1) ) ; 
             b=b+1 ;
        end
        %     + squeeze( R(idx1,idx2+2,:) ) + squeeze( R(idx1+1,idx2+3,:) )
        %     + squeeze( R(idx2+2,idx1,:) ) + squeeze( R(idx2+3,idx1+1,:) )

                            % subplot(2,2,3); 
        %                 plot(kc, abs(squeeze(C(idx1,idx2,:))) / temp_d);
        %                 title(sprintf('Periodic cross-correlation X(:,%d) and X(:,%d)', idx1, idx2)); 
        %                 xlabel('k'); ylabel(temp_str); grid on;
        % subplot(2,2,4); 
        %                 plot(kc, abs(squeeze(C(idx2,idx1,:))) / temp_d);
        %                 title(sprintf('Periodic cross-correlation X(:,%d) and X(:,%d)', idx2, idx1));
        %                 xlabel('k'); ylabel(temp_str); grid on;
        % 
        % %% psl and isl against sequence length
        % params    = [];
        % params.M  = M;
        % params.R  = r;
        %     
        % Apsl      = zeros(1, N_max); 
        % Aisl      = zeros(1, N_max);
        % Ppsl      = zeros(1, N_max);
        % Pisl      = zeros(1, N_max);
        % welchisl  = zeros(1, N_max);
        % welchApsl = zeros(1, N_max);
        % welchPpsl = zeros(1, N_max);
        % 
        % warning('Welch bound is not valid for M = 1 ');
        % for n = 2:N_max
        %     params.N     = n;
        %     X            = genSignal(ch_X, params);
        %     R            = correlation(X, 'a');
        %     C            = correlation(X, 'p');
        %     Apsl(n)      = psl(R); 
        %     Aisl(n)      = isl(R);
        %     Ppsl(n)      = psl(C); 
        %     Pisl(n)      = isl(C);
        %     sigma        = norm(X,'fro')^2/M;
        %     welchApsl(n) = sigma * sqrt((M-1)/(2*M*n-M-1));
        %     welchPpsl(n) = sigma * sqrt((M-1)/(M*n-1));
        %     welchisl(n)  = sigma^2 * M *(M-1);
        % end
        % 
        % % plots
        % figure(2); 
        % subplot(2,2,1); 
        %                 plot(Apsl);  hold on;
        %                 plot(welchApsl);  
        %                 title('Aperiodic PSL growth'); legend('PSL', 'Welch bound');
        %                 xlabel('Sequence length'); grid on; hold off;
        % subplot(2,2,2); 
        %                 plot(Aisl); hold on;
        %                 plot(welchisl);
        %                 title('Aperiodic ISL growth'); legend('ISL', 'Welch bound');
        %                 xlabel('Sequence length'); grid on;  hold off;
        % subplot(2,2,3); 
        %                 plot(Ppsl); hold on;
        %                 plot(welchPpsl);
        %                 title('Periodic PSL growth'); legend('PSL', 'Welch bound');
        %                 xlabel('Sequence length'); grid on;  hold off;
        % subplot(2,2,4); 
        %                 plot(Pisl); hold on;
        %                 plot(welchisl);
        %                 title('Periodic ISL growth'); legend('ISL', 'Welch bound');
        %                 xlabel('Sequence length'); grid on;  hold off;
        % 
        %                 
        % %% frequency and phase distribution of a signle sequence
        % params    = [];
        % params.M  = 1;
        % params.N  = N;
        % params.R  = r;
        % 
        % x         = genSignal(ch_x, params);
        % X         = fftshift(fft(x,N));
        % 
        % % plots
        % figure(3);
        % subplot(2,1,1); plot(20*log10(abs(X))); title('Magnitude of the frequency spectrum');
        %                 xlabel('Frequency (N)'); ylabel('magnitude (dB)'); grid on;
        % subplot(2,1,2); plot(angle(X)/pi); title('Phase of the frequency spectrum');
        %                 xlabel('Frequency (N)'); ylabel('phase / \pi'); grid on;            