function [R_t] = LW_PRE2012_RHS_S123_Jul2025(R_in_orig,N,Bo,a,la,da_fact);
        
    maxk=floor(N/2/da_fact);
    kgone=[maxk+1:N-maxk+1];
    zeropad=zeros(1,length(kgone));
    wavenumber_vector=[0:N/2 -N/2+1:-1]/la;
    wavenumber_vector(kgone)=zeropad;

    % Calculate fft(R)
    Rhat=fft(R_in_orig);
    Rhat(kgone)=zeropad;
    R=ifft(Rhat);
    
    % Calculate derivatives of R and Gamma
    Rhatz=Rhat.*(1i*wavenumber_vector);
    % Rhatzz=Rhat.*(1i*wavenumber_vector).^2;
    Rhatzzz=Rhat.*(1i*wavenumber_vector).^3;
    Rz=ifft(Rhatz);
    Rzzz=ifft(Rhatzzz);

    % Calculate needed powers (and inverses) of R
    Rhat2=fft(R.*R);
    Rhat2(kgone)=zeropad;
    R2=ifft(Rhat2);
    
    Rhat4=fft(R2.*R2);
    Rhat4(kgone)=zeropad;
    R4=ifft(Rhat4);
    
    Rinvhat=fft(1./R);
    Rinvhat(kgone)=zeropad;
    Rinv=ifft(Rinvhat);
    
    Rinv2hat=fft(Rinv.*Rinv);
    Rinv2hat(kgone)=zeropad;
    Rinv2=ifft(Rinv2hat);
    
    Rinv2Rzhat=fft(Rinv2.*Rz);
    Rinv2Rzhat(kgone)=zeropad;
    Rinv2Rz=ifft(Rinv2Rzhat);
    
    % Rinv4hat=fft(Rinv2.*Rinv2);
    % Rinv4hat(kgone)=zeropad;
    % Rinv4=ifft(Rinv4hat);
    
    % Rinv6hat=fft(Rinv2.*Rinv4);
    % Rinv6hat(kgone)=zeropad;
    % Rinv6=ifft(Rinv6hat);
    
    % Calculate log and some products
    logRhat=fft(log(R/a));
    logRhat(kgone)=zeropad;
    logR=ifft(logRhat);
    
    % R2logRhat=fft(R2.*logR);
    % R2logRhat(kgone)=zeropad;
    % R2logR=ifft(R2logRhat);
    
    R4logRhat=fft(R4.*logR);
    R4logRhat(kgone)=zeropad;
    R4logR=ifft(R4logRhat);
    
    % R2Rzzzhat=fft(R2.*Rzzz);
    % R2Rzzzhat(kgone)=zeropad;
    % R2Rzzz=ifft(R2Rzzzhat);
    
    % Calculate needed f's, P's, and wl (see report)
    f1=a^4-4*a^2.*R2+3*R4-4*R4logR;
    f1hat=fft(f1);
    f1hat(kgone)=zeropad;
    f1=ifft(f1hat);
    
    % f2=R2-a^2-2*R2logR-2*Ltilde*a^2+2*Ltilde.*R2;
    % % f2 is f2 from SOT23 and is 2*f1fromCOO14
    % f2hat=fft(f2);
    % f2hat(kgone)=zeropad;
    % f2=ifft(f2hat);
    % 
    % f3=a^4.*Rinv2-4*a^2+3*R2-4*R2logR+4*a^4*Ltilde.*Rinv2-8*a^2*Ltilde+4*Ltilde.*R2;
    % % f3 is f3 from SOT23 and is -f2fromCOO14
    % f3hat=fft(f3);
    % f3hat(kgone)=zeropad;
    % f3=ifft(f3hat);

    P1=f1.*(1-a^2/Bo*(Rinv2Rz+Rzzz));
    P1hat=fft(P1);
    P1hat(kgone)=zeropad;
    P1zhat=P1hat.*(1i*wavenumber_vector);
    P1z=ifft(P1zhat);

    % Calculate R_t and Gamma_t
    R_t = 1/16.*Rinv.*P1z;
    
    R_t_hat=fft(R_t);
    R_t_hat(kgone)=zeropad;
    R_t=ifft(R_t_hat);

end
