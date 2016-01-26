function ax=nndetector_vis_test(NTARGETS,TESTOUT,SPECTROGRAMS,TIMES,...
  TIME_WIN,TIME_WIN_STEPS,TRIGGER_THRESHOLDS,NTRAIN,NTEST,TIMESTEP,RNDSONGS,NMATCH)

SHOW_THRESHOLDS=true;
SORT_BY_ALIGNMENT=true;

disp('Creating spectral power image...');

% Create an image on which to superimpose the results...
power_img = squeeze((sum(SPECTROGRAMS, 2)));
power_img(find(isinf(power_img))) = 0;
power_img = power_img(RNDSONGS,:);
power_img = repmat(power_img / max(max(power_img)), [1 1 3]);

nsongs=size(SPECTROGRAMS,1);
AX=[];

for i = 1:NTARGETS
        AX(i)=subplot(NTARGETS+1,1,i+1);
        foo = reshape(TESTOUT(i,:,:), [], nsongs);
        barrr = zeros(TIME_WIN_STEPS-1, nsongs);

        if SHOW_THRESHOLDS
                img = power_img * 0.8;
                fooo = nndetector_trigger(foo', TRIGGER_THRESHOLDS(i), 0.1, TIMESTEP);
                fooo = [barrr' fooo];
                [val pos] = max(fooo,[],2);

                img(1:NTRAIN, :, 1) = img(1:NTRAIN, :, 1) - fooo(1:NTRAIN,:);
                img(1:NTRAIN, :, 2) = img(1:NTRAIN, :, 2) + fooo(1:NTRAIN,:);
                img(1:NTRAIN, :, 3) = img(1:NTRAIN, :, 3) + fooo(1:NTRAIN,:);
                img(NTRAIN+1:end, :, 1) = img(NTRAIN+1:end, :, 1) + fooo(NTRAIN+1:end,:);
                img(NTRAIN+1:end, :, 2) = img(NTRAIN+1:end, :, 2) - fooo(NTRAIN+1:end,:);
                img(NTRAIN+1:end, :, 3) = img(NTRAIN+1:end, :, 3) - fooo(NTRAIN+1:end,:);

                % color training, testing, and all negative examples

                img(1:NTRAIN, 1:TIME_WIN_STEPS, 1) = 0;
                img(1:NTRAIN, 1:TIME_WIN_STEPS, 2) = 1;
                img(1:NTRAIN, 1:TIME_WIN_STEPS, 3) = 1;

                img(NTRAIN+1:end, 1:TIME_WIN_STEPS, 1) = 1;
                img(NTRAIN+1:end, 1:TIME_WIN_STEPS, 2) = 0;
                img(NTRAIN+1:end, 1:TIME_WIN_STEPS, 3) = 0;

                %img(NMATCH+1:end, TIME_WIN_STEPS+1:TIME_WIN_STEPS*2, 1) = 1;
                %img(NMATCH+1:end, TIME_WIN_STEPS+1:TIME_WIN_STEPS*2, 2) = 0;
                %img(NMATCH+1:end, TIME_WIN_STEPS+1:TIME_WIN_STEPS*2, 3) = 1;

                if SORT_BY_ALIGNMENT
                        %[~, new_world_order] = sort(target_offsets);
                        [~, new_world_order] = sort(pos);
                        img = img(new_world_order,:,:);
                end

                image([TIMES(1) TIMES(end)], [1 nsongs], img);
        else
                barrr(:, 1:NTRAIN) = max(max(foo))/2;
                barrr(:, NTRAIN+1:end) = 3*max(max(foo))/4;
                foo = [barrr' foo'];
                imagesc([TIMES(1) TIMES(end)], [1 nsongs], foo);
        end
        xlabel('Time (ms)');
        ylabel('Song (random order)');
        if ~SORT_BY_ALIGNMENT
                text(time_window/2, NTRAIN/2, 'train', ...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Rotation', 90);
                text(time_window/2, NTRAIN+NTEST/2, 'test', ...
                        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Rotation', 90);
        end
end
