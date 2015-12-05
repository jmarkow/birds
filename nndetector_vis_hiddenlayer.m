function nndetector_vis_hiddenlayer(NET,FFT_TIME_SHIFT,TIME_WIN_STEPS,FREQ_RANGE,FREQ_RANGE_DS)
%
%
%

per_row=5;

if NET.numLayers > 1
  nunits=size(NET.IW{1},1);
  nrows=ceil(nunits/5);
  per_row=min(per_row,nunits);

  for i = 1:nunits
    subplot(nrows,per_row, i)
    imagesc([-TIME_WIN_STEPS:0]*FFT_TIME_SHIFT*1000, linspace(FREQ_RANGE(1), FREQ_RANGE(2), length(FREQ_RANGE_DS))/1000, ...
    reshape(NET.IW{1}(i,:), length(FREQ_RANGE_DS), TIME_WIN_STEPS));
    axis xy;
    axis off;
    if i == 1
      title('Hidden units');
    end

    if i==nunits
      hold on;
      xlimits=xlim();
      ylimits=ylim();
      left_edge=xlimits(2)+diff(xlimits)*.1;
      bot_edge=ylimits(1);
      h=line([left_edge left_edge],[bot_edge bot_edge+1]);
      h2=line([left_edge left_edge+10],[bot_edge bot_edge]);
      set(h,'clipping','off');
      set(h2,'clipping','off');
    end

  end
end
