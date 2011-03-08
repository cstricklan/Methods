function y = expandB(b, m)
tic
y=reshape(repmat(b',1,m)', length(b(:,1)), m*length(b(1,: )));
disp('using reshape/repmat')
toc
tic
y = kron(b, ones(1,m));
disp('using kron')
toc