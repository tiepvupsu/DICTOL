function [y,nz] = ompdenoise1(params,msgdelta)
%OMPDENOISE1 OMP denoising of 1-D signals.
%  OMPDENOISE1 denoises a 1-dimensional signal using OMP denoising. The
%  function syntax is identical to OMPDENOISE, but it runs significantly
%  faster on 1-D signals. OMPDENOISE1 requires somewhat more memory than
%  OMPDENOISE (approximately the size of the input signal), so if memory is
%  limited, OMPDENOISE can be used instead.
%
%  See also OMPDENOISE.


%  Ron Rubinstein
%  Computer Science Department
%  Technion, Haifa 32000 Israel
%  ronrubin@cs
%
%  August 2009


% parse input arguments %

x = params.x(:);
D = params.dict;
blocksize = params.blocksize;


% maxval %
if (isfield(params,'maxval'))
  maxval = params.maxval;
else
  maxval = 1;
end


% gain %
if (isfield(params,'gain'))
  gain = params.gain;
else
  gain = 1.15;
end


% maxatoms %
if (isfield(params,'maxatoms'))
  maxatoms = params.maxatoms;
else
  maxatoms = floor(blocksize/2);
end


% stepsize %
if (isfield(params,'stepsize'))
  stepsize = params.stepsize;
else
  stepsize = 1;
end
if (any(stepsize<1))
  error('Invalid step size.');
end


% noise mode %
if (isfield(params,'noisemode'))
  switch lower(params.noisemode)
    case 'psnr'
      sigma = maxval / 10^(params.psnr/20);
    case 'sigma'
      sigma = params.sigma;
    otherwise
      error('Invalid noise mode specified');
  end
elseif (isfield(params,'sigma'))
  sigma = params.sigma;
elseif (isfield(params,'psnr'))
  sigma = maxval / 10^(params.psnr/20);
else
  error('Noise strength not specified');
end


% lambda %
if (isfield(params,'lambda'))
  lambda = params.lambda;
else
  lambda = maxval/(10*sigma);
end


% msgdelta %
if (nargin <2)
  msgdelta = 5;
end
if (msgdelta<=0)
  msgdelta = -1;
end


epsilon = sqrt(blocksize) * sigma * gain;   % target error for omp


MEM_LOW = 1;
MEM_NORMAL = 2;
MEM_HIGH = 3;

if (isfield(params,'memusage'))
  switch lower(params.memusage)
    case 'low'
      memusage = MEM_LOW;
    case 'normal'
      memusage = MEM_NORMAL;
    case 'high'
      memusage = MEM_HIGH;
    otherwise
      error('Invalid memory usage mode');
  end
else
  memusage = MEM_NORMAL;
end


% compute G %

G = [];
if (memusage >= MEM_NORMAL)
  G = D'*D;
end


% verify dictionary normalization %

if (isempty(G))
  atomnorms = sum(D.*D);
else
  atomnorms = diag(G);
end
if (any(abs(atomnorms-1) > 1e-2))
  error('Dictionary columns must be normalized to unit length');
end


% denoise the signal %


% process the signal in batches to conserve memory
% choose batchsize so im2colstep returns a matrix of approximately the same
% size as the signal
batchsize = ceil(length(x)*stepsize/blocksize + blocksize);

y = zeros(size(x));
ids = 1:min(batchsize,length(x));
nz = 0;

blocknum = floor((length(x)-blocksize)/stepsize) + 1;
processedblocks = 0;
tid = timerinit('ompdenoise', blocknum);
while (length(ids)>=blocksize)

  % extract the signal blocks
  blocks = im2colstep(x(ids),[blocksize 1],[stepsize 1]);

  % remove DC
  [blocks, dc] = remove_dc(blocks,'columns');

  % denoise the blocks
  if (memusage == MEM_LOW)
    gamma = omp2(D,blocks,[],epsilon,'maxatoms',maxatoms,'checkdict','off');
  else
    gamma = omp2(D'*blocks,sum(blocks.*blocks),G,epsilon,'maxatoms',maxatoms,'checkdict','off');
  end
  nz = nz + nnz(gamma);
  cleanblocks = add_dc(D*gamma, dc, 'columns');

  y(ids) = y(ids) + col2imstep(cleanblocks, [length(ids) 1], [blocksize 1], [stepsize 1]);
  ids = ids + floor((batchsize-blocksize)/stepsize)*stepsize + stepsize;
  if (ids(end)>length(x))
    ids = ids(ids<=length(x));
  end
  
  % display status
  if (msgdelta>0)
    processedblocks = processedblocks + size(blocks,2);
    timereta(tid, processedblocks, msgdelta);
  end

end

if (msgdelta>0)
  timereta(tid, blocknum);
end


cnt = countcover(size(x),[blocksize 1],[stepsize 1]);
y = (y+lambda*x)./(cnt + lambda);
y = reshape(y,size(params.x));
