function ParForUpdateWaitbar(minmax)

global loopcount
global totloopcount
global strCR

if nargin<1||isempty(minmax)
    minmax = [0 100];
end

strPercentageLength = 10;
strDotsMaximum = 10;
 
pctg = (minmax(2)-minmax(1))*loopcount/totloopcount+minmax(1);
pctg = floor(pctg);
percentageOut = [num2str(pctg) '%%'];
percentageOut = [percentageOut repmat(' ',1,strPercentageLength-length(percentageOut)-1)];
nDots = floor(pctg/100*strDotsMaximum);
dotOut = ['[' repmat('.',1,nDots) repmat(' ',1,strDotsMaximum-nDots) ']'];
strOut = [percentageOut dotOut];

% Print it on the screen
if strCR == -1,
    % Don't do carriage return during first run
    fprintf(strOut);
else
    % Do it during all the other runs
    fprintf([strCR strOut]);
end

% Update carriage return
strCR = repmat('\b',1,length(strOut)-1);

loopcount = loopcount + 1;
end