function nndetector_vis_hiddenlayer(NET)
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
    imagesc([-NET.userdata.time_window_steps:0]*NET.userdata.fft_time_shift/NET.userdata.samplerate*1000, ...
      linspace(NET.userdata.freq_range(1), NET.userdata.freq_range(2), length(NET.userdata.freq_range))/1e3, ...
      reshape(NET.IW{1}(i,:), length(NET.userdata.freq_range_ds), NET.userdata.time_window_steps));
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
