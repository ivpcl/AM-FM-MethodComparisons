% -------------------------------------------------------------------------
% AM-FM Methods Experimental Comparisons
% -------------------------------------------------------------------------
% Authors: Kyriacos Constantinou, Ioannis Constantinou,
% Marios S. Pattichis and Constantinos S. Pattichis
% csp7ck1@cs.ucy.ac.cy, ioannis@istognossis.com
% pattichi@unm.edu, pattichi@cs.ucy.ac.cy
%
% This file is available at:
% https://ivpcl.unm.edu/ivpclpages/Research/amfm/amfm.php
% http://www.ehealthlab.cs.ucy.ac.cy/index.php/facilities/32-software/225-toolboxes
%
% This code provide comparisons between the following AM-FM methods:
% (i) Multiscale AM-FM method using a Gabor Filterbank,
% (ii) Multidimensional Ensemble Empirical Mode Decomposition method, and
% (iii) Variational Mode Decomposition for the modified Lagrangian.
%
% Initial release 07 November 2019 (c) 2019
%
% When using this code, please give credit to the following paper:
% -------------------------------------------------------------------------
% K. Constantinou, I. Constantinou, M. S. Pattichis and C. S. Pattichis,
% Medical Image Analysis Using AM-FM Models and Methods,
% IEEE Reviews on Biomedical Engineering, revised version submitted Nov. 2019.
%
% Note: This code generates the analysis for Fig. 6. 
% A comparison of AM-FM component estimation for three different methods 
% (of the aforementioned paper).
% -------------------------------------------------------------------------
% This version was run on Matlab version 9.5.0 (R2018b)
% -------------------------------------------------------------------------
%% Load synthetic image
% -------------------------------------------------------------------------
img = syntheticImage();
% -------------------------------------------------------------------------
%%  Multiscale AM-FM using Gabor Filerbank
% -------------------------------------------------------------------------
tic
filters = setFilterbank();
hImg    = hilbert(img);
AMFM    = {};
for k = 1:length(filters)
    filter          = filters{k};
    filterImg       = filter2(filter, hImg);
    [IA,IP,IFx,IFy] = calculateAMFM(filterImg);
    AMFMRes         = struct('IA',IA , 'IP', IP, 'IFx', IFx, 'IFy', IFy);
    AMFM            = cat(1, AMFM, AMFMRes);
end

% Multiscale Dominant Component Analysis
high = {}; med = {}; low = {}; dc = {};
for i=1:length(AMFM)
    switch i
        case num2cell(1:8)
            high = cat(1,high, AMFM(i));
        case num2cell(8:24)
            med = cat(1,med, AMFM(i));
        case num2cell(25:40)
            low = cat(1,low, AMFM(i));
        case num2cell(41)
            dc = cat (1,dc, AMFM(i));
    end
end

[ IAl, IPl, IFxl, IFyl ] = dca (low);
% Thresholding IA
IAl = (IAl > prctile(IAl(:),50)).*IAl;
reconstructionImgDCAl = real(squeeze(IAl(:,:).*cos(IPl(:,:))));


[ IAm, IPm, IFxm, IFym ] = dca (med);
% Thresholding IA
IAm = (IAm > prctile(IAm(:),96)).*IAm;
reconstructionImgDCAm = real(squeeze(IAm(:,:).*cos(IPm(:,:))));


[ IAh, IPh, IFxh, IFyh ] = dca (high);
% Thresholding IA
IAh = (IAh > prctile(IAh(:),96)).*IAh;
reconstructionImgDCAh = real(squeeze(IAh(:,:).*cos(IPh(:,:))));

[ IAdc, IPdc, IFxdc, IFydc ] = dca(dc);
reconstructionImgDCAdc = real(squeeze(IAdc(:,:).*cos(IPdc(:,:))));
toc

figure, imagesc(reconstructionImgDCAl), colormap gray, title ('DCA Low Scale Reconstruction'),axis equal off;
figure, imagesc(reconstructionImgDCAm), colormap gray, title ('DCA Medium Scale Reconstruction'),axis equal off;
figure, imagesc(reconstructionImgDCAh), colormap gray, title ('DCA High Scale Reconstruction'),axis equal off;
figure, imagesc(reconstructionImgDCAdc), colormap gray, title ('DCA DC Reconstruction'),axis equal off;

function filters = setFilterbank()
%   Filterbank Parameters
lambda0      = 2;
orientations = 8;
scales       = 5;
gamma        = 1;
phase        = 0;
bandwidth    = 1;
overlapIndex = 0.5;
offset       = 0;
theta        = offset:(pi/orientations):(pi - pi/orientations + offset);
lambda       = lambda0;
filters      = {};

for sc_index = scales:-1:1
    lambda0 = lambda;
    for th=1:length(theta)
        [result,sig] = gaborKernel2D(theta(th), lambda, gamma,...
            bandwidth, phase, 1/overlapIndex);
        filters      = cat(1, filters, result);
    end
    lambda = lambda0*2^bandwidth;
end
%   Add DC Filter
f1      = 2*pi/lambda;
result  = gaussianFunction(f1, sig, overlapIndex );
filters = cat(1, filters, result);

    function [result, sigma] = gaborKernel2D(theta, lambda, gamma,...
            bandwidth, phase, overlapIndex)
        qFactor   = (1/pi) * sqrt( (log(overlapIndex)/2) ) * ...
            ( (2^bandwidth + 1) / (2^bandwidth - 1) );
        sigma     = lambda*qFactor;
        n         = ceil(4*sigma);
        [x,y]     = meshgrid(-n:n+1,-n:n+1);
        xTheta    = x * cos(theta) + y * sin(theta);
        yTheta    = -x * sin(theta) + y * cos(theta);
        gaussian  = exp(-(( xTheta.^2) + gamma^2.* (yTheta.^2))./(2*sigma^2));
        res       = gaussian.* cos(2*pi*xTheta/lambda +  phase);
        maxFft    = max(max(abs(fft2(res(:,:)))));
        normalize = fft2(res(:,:))./maxFft;
        result    = real(ifft2(normalize));
    end

    function res = gaussianFunction(f0, s0, overlapIndex)
        over  = sqrt(2*log(1/overlapIndex));
        sigma = s0*over/(s0*f0 - over);
        n     = ceil(2*sigma);
        [x,y] = meshgrid(-n:n+1,-n:n+1);
        res   = exp(-1/2*(x.^2 + y.^2)./sigma^2);
        res   = res ./ (sum(res(:)));
    end
end

function [IA, IP, IFx, IFy] = calculateAMFM(filterImg)
[M,N]  = size(filterImg);
% Calculate IA
IA     = abs(filterImg);
% Calculate Phase
IP     = angle(filterImg);
% IA normilized
IANorm = filterImg ./IA;
IFx    = zeros(M,N);
IFy    = zeros(M,N);
for i=2:M-1
    for j=2:N-1
        IFx(i,j) = abs(acos(real((IANorm(i+1,j) + IANorm(i-1,j)) ./ (2*IANorm(i,j)))));
        IFy(i,j) = abs(acos(real((IANorm(i,j+1) + IANorm(i,j-1)) ./ (2*IANorm(i,j)))));
    end
end
end

function [ IA, IP, IFx, IFy ] = dca(band)
IA    = [];
IP    = [];
[w,h] = size(band{1}.IA);
for i=1:w
    for j=1:h
        % find the max for (i,j) position
        pos = 1;
        temp = band{pos}.IA(i,j);
        for l=1:length(band)
            if temp < band{l}.IA(i,j)
                pos = l;
                temp = band{l}.IA(i,j);
            end
        end
        IA(i,j)  = temp;
        IP(i,j)  = band{pos}.IP(i,j);
        IFx(i,j) = band{pos}.IFx(i,j);
        IFy(i,j) = band{pos}.IFy(i,j);
    end
end
end

% -------------------------------------------------------------------------
%% Multidimensional Ensamble Empirical Mode Decomposition
%
% Matlab code for the method taken from the Apendix of
%
% Z. WU, N. E. HUANG, and X. CHEN, The Multi-Dimensional Ensemble
% Empirical Mode Decomposition Method, Adv. Adapt. Data Anal., vol. 01,
% no. 03, pp. 339?372, 2009.
%
% -------------------------------------------------------------------------
% Example initial parameters as follow
% -------------------------------------------------------------------------
%
% number of modes 3
%
% -------------------------------------------------------------------------
%% Variational Mode Decomposition
%
% Matlab code taken from
%
% http://www.math.montana.edu/dzosso/code/
%
% -------------------------------------------------------------------------
% References
% -------------------------------------------------------------------------
%
% [1] D. Zosso, K. Dragomiretskiy, A.L. Bertozzi, P.S. Weiss, Two-Dimensional
% Compact Variational Mode Decomposition, Journal of Mathematical Imaging
% and Vision, 58(2):294–320, 2017.
% DOI:10.1007/s10851-017-0710-z
%
% [2] K. Dragomiretskiy, D. Zosso, Variational Mode Decomposition, IEEE Trans.
% on Signal Processing, 62(3):531-544, 2014.
% DOI:10.1109/TSP.2013.2288675
%
% [3] K. Dragomiretskiy, D. Zosso, Two-Dimensional Variational Mode
% Decomposition, EMMCVPR 2015, Hong Kong, LNCS 8932:197-208, 2015.
% DOI:10.1007/978-3-319-14612-6_15
%
% -------------------------------------------------------------------------
% Example initial parameters as follow
% -------------------------------------------------------------------------
%
% alpha = 1000;       % bandwidth constraint
% beta = 0.5;         % L1 penalty on A (compact support)
% gamma = 500;         % TV weight
% delta = inf;        % threshold value for artifact detection
% rho = 10;           % fidelity weight
% rho_k = 10;        % u/v coupling weight
% tau = 2.5;          % gradient ascent step size for fid. lagrange multipliers
% tau_k = 2.5;        % gradient ascent step size for u/v lagrange multipliers
% t = 1.5;             % gradient descent step size for A heat propogation
%
% K = 4;              % number of modes
% M = 1;              % number of submodes (per mode)
% DC = 1;             % includes DC part (first mode at DC)
% init = 0;           % initialize omegas
%
% u_tol = 1e-10;      % tolerance for u convergence
% A_tol = 1e-4;    % tolerance for A convergence
% omega_tol = 1e-10;  % tolerance for omega convergence
%
% N = 130;             % max number of iterations
%
% A_phase = [100,Inf]; % A propogation phases [a,b]
%
% signal = img;
%
% -------------------------------------------------------------------------

function img = syntheticImage()
% Initialization parameters
N = 128; M = N;
a = N/4; b = M/4;
x0                = [5*N/8 5*N/8 N/4];
y0                = [3*M/8 3*M/4 M/2];
theta             = [3*pi/4 3*pi/4 pi/4];
rho               = [pi/16 pi/2 3*pi/4];
sigma             = [8/(9*pi) 1/(2*pi) 1/(2*pi)];
img               = double(zeros(N,M));
for i=1:length(x0)
    [omega1, omega2] = pol2cart(theta(i), rho(i));
    comp = buildAMFM(N, M, omega1, omega2, sigma(i), x0(i), y0(i), a, b);
    figure; imagesc(comp), colormap gray, axis off image;
    img = img + comp;
end
noise = randn(size(img))*5;
figure;imagesc(noise), colormap gray, axis off image;
img = img + noise;
figure; imagesc(img), colormap gray, axis off image;

    function comp = buildAMFM(N, M, omega1, omega2, sigma, x0, y0, a, b)
        [origX, origY] = meshgrid(0:1:N-1, 0:1:M-1);
        x = (origX - x0) ./ (N/2);
        y = (origY - y0) ./ (M/2);
        % Build AM
        amp = 100*exp(-1/2*( x.^2 + y.^2)./sigma^2);
        % Build FM
        phi = (a)*(x).^2 + (omega1)*(origX - x0) + (b)*(y).^2 + (omega2)*(origY - y0);
        % AM-FM component
        comp = (amp).*(cos(phi));
    end
end