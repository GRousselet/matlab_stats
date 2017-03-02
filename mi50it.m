function it = mi50it(data,xf,pts,plotit)
% it = mi50it(data,xf,pts,plotit)
%
% Computes the 50% integration time of mutual information data.
% 50it is a measure of processing speed that alleviates the need to consider peaks,
% by considering the shape of a time-course in a time-window of interest.
% The measure was introduced in this paper:
% Rousselet, G.A., Gaspar, C.M., Pernet, C.R., Husk, J.S., Bennett, P.J. & Sekuler, A.B. (2010)
% Healthy aging delays scalp EEG sensitivity to noise in a face discrimination task.
% Front Psychol, 1, 19.
% http://www.frontiersin.org/perception_science/10.3389/fpsyg.2010.00019/abstract
% The current version is slightly different, starting with a baseline correction
% so as to not integrate MI values around baseline level.
%
% data      = a column vector or time x participant matrix of positive only observations, e.g. mutual information
% xf        = time vector - default -300:2:600
% pts       = points at which to estimate the time-course - default [0 500]
% plotit    = 1 to get a figure, 0 otherwise

% Copyright (C) 2016 Guillaume Rousselet - University of Glasgow

if nargin < 2
    xf = -300:2:600;
end
if nargin < 3
    pts = [0 500];
end
if nargin < 4
    plotit = 0;
end

[nf,np] = size(data);
frame_selec = xf>=pts(1) & xf<=pts(2);
xf50it = xf(frame_selec);
nf50it = numel(xf50it);

% baseline correction
% data = data - repmat( median(data(xf<0,:),1), [nf 1] );

% cumulative sum
data = data(frame_selec,:);
data = cumsum( data, 1 );

% normalization
data = data - repmat(min(data,[],1), [nf50it 1]); % figure;plot(Xf0500,cs1)
data = data ./ repmat(max(data,[],1), [nf50it 1]); % figure;plot(Xf0500,cs1)
it = zeros(np,1);
for P = 1:np % for every participant
    it(P) = interp1(data(:,P),xf50it,.5,'linear');
end

if plotit == 1
    
    cc = parula(np);
    figure('Color','w','NumberTitle','off','Name','50% integration time')
    hold on
    for P=1:np
        plot(xf50it,data(:,P),'Color',cc(P,:),'LineWidth',2)
        plot([it(P) it(P)],[0 1],'Color',cc(P,:),'LineWidth',1)
    end
    set(gca,'LineWidth',2,'YLim',[0 1],'XLim',pts,'XTick',pts(1):50:pts(2))
    set(gca,'FontSize',14,'Layer','Top')
    xlabel('Time in ms','FontSize',16')
    ylabel('Normalised MI','FontSize',16')
    box on
    
end
