%open octave and add mvapack
addpath('/home/louic/Desktop/mvapack');

%icoshift data make sure to make ppm a separate file

spectra  = load('/home/louic/Desktop/mouse_pellets/150716_media/mario_media/spect_raw');
ppm = load('/home/louic/Desktop/mouse_pellets/150716_media/mario_media/ppm');
plot(ppm, spectra, 'b');
hold on;
segments = [0.85 1.15; 1.22 1.25; 1.25 1.40; 1.40 1.45; 1.45 1.50;...
	    1.55 1.60; 1.60 1.77; 1.77 1.80; 1.86 1.97; 1.97 2.10;...
	    2.10 2.20; 2.33 2.38; 2.38 2.42; 2.42 2.47; 2.47 2.57;...
	    2.62 2.71; 2.71 2.79; 2.79 2.81; 2.81 2.90; 2.91 2.98;...
	    2.98 3.08; 3.08 3.20; 3.20 3.31; 3.31 3.45; 3.45 3.63;...
	    3.63 3.70; 3.70 3.82; 3.82 3.87; 3.87 3.92; 3.92 3.97;...
	    3.97 4.03; 4.03 4.08; 4.08 4.16; 4.16 4.23; 4.23 4.29;...
	    4.29 4.38; 4.38 4.43; 4.43 4.60; 4.63 4.70; 5.00 5.10;...
	    5.20 5.25; 5.29 5.31; 5.31 5.35; 5.35 5.38; 5.38 5.42;...
	    5.45 5.48; 5.48 5.60; 6.85 6.95; 7.05 7.10; 7.15 7.22;...
	    7.30 7.50; 7.75 7.79; 7.79 7.85; 8.45 8.50];

%exectution 
%global then segmented
%if you get an error it could be because there are spectra that go below -1
shifted = icoshift(spectra, ppm, segments, cofirst =true); %cofirst will perform a global alignment before icoshifting
%X.s.data = icoshift(X.data, X.s.ppm, cofirst = true);

save "/home/louic/Desktop/mouse_pellets/150716_media/mario_media/spect_icoshift" shifted;
%segmented
%shifted = icoshift(spectra, ppm, segments);
plot(ppm, shifted, 'r');
hold off;
%save
save "/home/louic/Desktop/media_nmr/nmr/screen/spect_noppm_icoshift" shifted;
save "/home/louic/Desktop/media_nmr/nmr/screen/spect_noppm_icoshift2" shifted2;
X.pca.data = shifted;
X.ppm = ppm;
% normalize the entire dataset for pca.
X.pca.data = pqnorm(X.pca.data);

% align the entire dataset for pca analysis.
X.pca.wbin = 0.003;
[X.pca.data, X.pca.ppm, X.pca.widths] = ...
  binadapt(X.pca.data, X.ppm, F.parms, X.pca.wbin);
%remove noise
X.pca.noise = [findnearest(X.pca.ppm, 5.20), findnearest(X.pca.ppm, 4.58)];
[X.pca.data, X.pca.ppm] = rmnoise(X.pca.data, X.pca.ppm, X.pca.noise);

% build a pca model.
mdl.pca = pca(X.pca.data);
mdl.pca = addclasses(mdl.pca, cls.pca.Y);
mdl.pca = addlabels(mdl.pca, cls.labels);

scoresplot(pcaMdl, 3);

% build an lda model of the pca scores.
mdl.lda = lda(scores(mdl.pca), cls.pca.Y);
mdl.lda = addlabels(mdl.lda, cls.labels);

%plot pca
