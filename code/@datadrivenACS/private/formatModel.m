function dd = formatModel(dd)

if isempty(dd.Model.Plant)
    error('Must specify a model')
end

G = dd.Model.Plant;

lenG= length(G);
G = reshape(G,[lenG,1]);

% check for sampling grid size
if (isa(G,'tf')) || (isa(G,'ss')) || (isa(G,'zpk') || (isa(G,'idpoly')) || (isa(G,'frd')) || (isa(G,'idfrd')))
    if ~isempty(fields(G.SamplingGrid))
        dd.M.ndim = numel(fieldnames(G.SamplingGrid));
    else
        G.SamplingGrid = struct('tmp',zeros(1,1,lenG));
        dd.M.ndim = 1;
    end
end

dd.M.nmod = length(G); % number of different models
dd.M.G = @(w,mod) squeeze(freqresp(G(:,:,mod),w));

cellSamplingGrid = struct2cell(G.SamplingGrid);
cellSamplingGrid = reshape(cellSamplingGrid,[1,numel(cellSamplingGrid)]);
dd.M.SamplingGrid = cell2mat(cellSamplingGrid); % SamplingGrid used for LPV


% Build frequency vector if not specified
if ~isempty(dd.Model.Frequency)
    % user specified
    dd.M.freqs = dd.Model.Frequency;
else
    if (isa(G,'frd') || isa(G,'idfrd'))
        % if frd, it G.Frequency (and trim if needed)
        Wtmp = G.Frequency;
        
        if G.Ts > 0
            dd.M.freqs = Wtmp(Wtmp<=pi/G.Ts);
        else
            dd.M.freqs = Wtmp;
        end
    else
        % retrieve frequency vector from Nyquist command, not advised
        [~,~,freqs] = nyquist(G(:,:,1));
        dd.M.freqs = freqs;
    end
end

%%