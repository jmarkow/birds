function [ stats ] = optimise_network_output_unit_trigger_thresholds(NET,NNSETX,NNSETY,WEIGHTS,SLOP,THRESHOLDS)
% slop to match (.005)?

if nargin<6
  THRESHOLDS=[];
end

if nargin<5
  SLOP=[];
end

if nargin<4
  WEIGHTS=[.5 .5];
end

t_per_step=NET.userdata.time_window/NET.userdata.time_window_steps;

if ~isempty(SLOP)
  slop_smps=round(SLOP/t_per_step)
  if slop_smps>1
    NNSETY=conv(NNSETY,ones(1,slop_smps),'same');
  end
end

activation=sim(NET,NNSETX);

if isempty(THRESHOLDS)
  THRESHOLDS=linspace(min(activation),max(activation),1e3);
end

% build roc curve

stats.accuracy=zeros(1,length(THRESHOLDS));
stats.fpr=zeros(size(stats.accuracy));
stats.tpr=zeros(size(stats.accuracy));
stats.tnr=zeros(size(stats.accuracy));
stats.f1=zeros(size(stats.accuracy));
stats.youden=zeros(size(stats.accuracy));
stats.thresholds=THRESHOLDS;
stats.activation=activation;

NNSETY=NNSETY>0;
condition_positive=sum(NNSETY)
condition_negative=sum(~NNSETY)
total=length(NNSETY)

for i=1:length(THRESHOLDS)

    prediction=activation>=THRESHOLDS(i);

    tp=sum(prediction&NNSETY);
    fp=sum(prediction&~NNSETY);
    tn=sum(~prediction&~NNSETY);
    fn=sum(~prediction&NNSETY);

    stats.accuracy(i)=(tp+tn)/total;
    stats.weighted_cost(i)=WEIGHTS(1)*fp+WEIGHTS(2)*fn;
    stats.tpr(i)=tp/(tp+fn);
    stats.fpr(i)=fp/(fp+tn);
    stats.tnr(i)=tn/(fp+tn);
    stats.f1(i)=2*tp/(2*tp+fp+fn);
    stats.youden(i)=(stats.tpr(i)+stats.tnr(i)-1);

end
