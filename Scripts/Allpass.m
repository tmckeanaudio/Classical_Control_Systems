%% Allpass System Lecture
% $$G(s) = \frac{1 - as}{1 + as}$$

s = tf('s'); a = 1/5;
G = (1 - a*s)/(1 + a*s)
figure, bode(G)
P = 1/(s + 1)
k = 1; C = k;
Y = P*G*k 
G = (-k/(s+1))*(s-1/a)/(s+1/a)
figure, bode(G), margin(G)
pzmap(G)
nyquist(G)
rlocus(G)
[k,poles] = rlocfind(G)
%% 
% Thus, k is about 6 so our range can only be within 
% 
% $$0<k<6$$