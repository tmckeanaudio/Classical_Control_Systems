%% Nyquist Stability Criterion
%% 
% * Derive Nyquist plots from Transfer Functions $G(s)$ where $s = \sigma + 
% j\omega$

s = tf('s');
G = 1/((s+2)*(s+4)*(s+6))
figure, nyquist(G), xlim([-0.005 .021])
syms s w
G = 1/((s+2)*(s+4)*(s+6))
G_w = subs(G,s,j*w)
W = -100:.01:100; 
Nyq = eval(subs(G_w,w,W));
x = real(Nyq);
y = imag(Nyq);
figure, plot(x,y)
%% 
% Now we can solve for the Gain margin: 
% 
% Want to find the points when our imaginary values change sign

tmp = sign(y);
d = diff(tmp);
indexes = find(d ~= 0)
points = W(indexes), a = points(4)
GainFactor = 1/abs(x(indexes(4)))
%% 
% Thus our Gain Factor for this transfer function is about 480
%% Assessing Gain and Phase Margins

f = imread('StabilityExample.png'); imshow(f)
G = tf([0.5 1.3],[1 1.2 1.6 0])
T = feedback(G,1)
p = pole(T)
p1mag = abs(p(1)), p1pha = angle(p(1))*(180/pi)
p2mag = abs(p(2)), p2pha = angle(p(2))*(180/pi)
p3mag = abs(p(3)), p3pha = angle(p(3))*(180/pi)
%% 
% Root Locus of G:

rlocus(G)
%% 
% Using the Root Locus, we want to assess how much the gain can change before 
% stability to the system is lost

[k,poles] = rlocfind(G)
%% 
% Clicking on the imaginary axis of the graph we see that when $k = 2.7$ is 
% when the system becomes unstable.
% 
% Thus our gain can be within the range of $0 < k < 2.7 $
% Gain and Phase Margins
% The phase margin measures how much phase variation is needed at the gain crossover 
% frequency to lose stability. Similarly, the gain margin measures what relative 
% gain variation is needed at the gain crossover frequency to lose stability. 
% Together, these two numbers give an estimate of the "safety margin" for closed-loop 
% stability. The smaller the stability margins, the more fragile stability is

bode(G), grid, margin(G)