function  [mean_u]= channel_model (TxLoc, RxLoc)

% Calculate  RSS mean value.

  d0 = 1;% Reference Distance
  n = 3; % n is the path loss exponent
  Gt = 1;
  Gr = 1;
  wavelen = (3*10^8)/(9*10^9); % lembda = c/f for 900 MHz band
  TxPow_L = 1000;   % equal power of 1W (1000mW)
  TxPow_dB = 10*log10(TxPow_L);
  PL_Ref_dB = -10*log10(Gt*Gr*wavelen^2./(4*pi*d0).^2);    % reference distance:1

    for i = 1 : size (RxLoc,1)
      d(i) = norm (TxLoc - RxLoc(i,:));
      PL_d_dB(i) = PL_Ref_dB + 10*n*log10(d(i)/d0);                  % no shadowing
      RSS_mean(i) = TxPow_dB - PL_d_dB(i);
    end
    
  mean_u=RSS_mean;

