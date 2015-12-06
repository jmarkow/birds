function nndetector_learn_export_wav(MIC_DATA,SAMPLE_RATE,TARGETS,SPREAD,FILENAME)
%
%
%
%
%
%

RIGHT_CHANNEL=zeros(size(MIC_DATA));

for i=1:length(TARGETS)
  RIGHT_CHANNEL(TARGETS(i),i)=1;
  RIGHT_CHANNEL(:,i)=conv(RIGHT_CHANNEL(:,i),ones(SPREAD,1),'same');
end

audiowrite([ FILENAME ],[MIC_DATA(:) RIGHT_CHANNEL(:)],SAMPLE_RATE);
