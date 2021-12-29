f = 200;
tt = [0:1/(20*f):1];
z = exp(i*2*pi*f*tt); 
subplot(211) 
plot(real(z))
title('Real part of exp(j*2*pi*200*tt)')
subplot(212)
plot(imag(z))
title('Imaginary part of exp(j*2*pi*200*tt)')