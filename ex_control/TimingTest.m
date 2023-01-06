function [miss] = TimingTest(allCodes,varargin)

p = inputParser;
p.addOptional('flag',false,@islogical);
p.parse(varargin{:});
flag = p.Results.flag;

miss = [];
for i = 1:length(allCodes)
    if ~isempty(find(allCodes{1,i}.codes(:,1)==251, 1))
        miss  = [miss i];
    end
end
fprintf('%d completed trials with %d frame drops\n',length(allCodes)-length(miss),length(miss))
fprintf('Percentage trials with frame drops = %f%%\n',length(miss)./length(allCodes)*100);
if flag
    figure;hold on;
    plot(miss,1,'r*')
    xlim([1 length(allCodes)]);
    text(miss,ones(size(miss)),num2cell(miss),'VerticalAlignment','bottom','HorizontalAlignment','right');
    title(sprintf('Dropped %d times in %d trials',length(miss),length(allCodes)-length(miss)));
end

