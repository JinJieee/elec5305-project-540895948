function study = placement_study(folder, source, responses, azimuth, t, route, fs, hop, N)
% Locations are defined relative to the listener at turn onset, not trial start.
bearings = -150:30:180;
distances = [30, 100, 300];
[position, heading] = listener_path(t, 1, route);
entry = [route.speed*(route.turnStart-route.walkStart), 0];
endTime = route.turnStart + route.turnDuration;
[~, first] = min(abs(t-route.turnStart));
[~, last] = min(abs(t-endTime));
values = zeros(numel(bearings)*numel(distances), 12);
worldCues = cell(numel(distances), numel(bearings));
headCues = worldCues;
angles = zeros(numel(t), numel(bearings), numel(distances));
row = 0;
for d = 1:numel(distances)
    for b = 1:numel(bearings)
        row = row + 1;
        fprintf('Placement %d/%d: %g deg, %g m\n', row, size(values,1), bearings(b), distances(d));
        alpha = bearings(b)*pi/180;
        location = entry + distances(d)*[cos(alpha), sin(alpha)];
        angle = unwrap(atan2(location(2)-position(:,2), location(1)-position(:,1)) - heading);
        assert(abs(atan2(sin(angle(first)-alpha),cos(angle(first)-alpha))) < 1e-10);
        assert(abs(norm(location-entry)-distances(d)) < 1e-10);
        world = render_binaural(source, angle, responses, azimuth, hop);
        head = render_binaural(source, alpha*ones(size(t)), responses, azimuth, hop);
        wc = measure_cues(world, fs, N);
        hc = measure_cues(head, fs, N);
        use = wc.Time_s >= 8.1 & wc.Time_s <= endTime-0.1;
        ildRange = max(wc.ILD_dB(use))-min(wc.ILD_dB(use));
        itdRange = max(wc.ITD_us(use))-min(wc.ITD_us(use));
        headILD = max(hc.ILD_dB(use))-min(hc.ILD_dB(use));
        headITD = max(hc.ITD_us(use))-min(hc.ITD_us(use));
        ildDifference = sqrt(mean((wc.ILD_dB(use)-hc.ILD_dB(use)).^2));
        itdDifference = sqrt(mean((wc.ITD_us(use)-hc.ITD_us(use)).^2));
        change = (angle(last)-angle(first))*180/pi;
        % Pure rotation gives -90 deg; the residual is the translation contribution.
        values(row,:) = [bearings(b), distances(d), location, change, change+90, ...
            ildRange, itdRange, headILD, headITD, ildDifference, itdDifference];
        worldCues{d,b} = wc;
        headCues{d,b} = hc;
        angles(:,b,d) = angle;
    end
end
summary = array2table(values, 'VariableNames', {'Bearing_deg','Distance_m','SourceX_m', ...
    'SourceY_m','AzimuthChange_deg','TranslationChange_deg','ILDRange_dB','ITDRange_us', ...
    'HeadILDRange_dB','HeadITDRange_us','PairedILDRMS_dB','PairedITDRMS_us'});
assert(all(isfinite(values(:))));
for b = 1:numel(bearings)
    rows = summary.Bearing_deg == bearings(b);
    assert(max(summary.HeadITDRange_us(rows))-min(summary.HeadITDRange_us(rows)) < 1e-9);
end
writetable(summary, fullfile(folder,'results','placement_study.csv'));
study = struct('summary',summary,'worldCues',{worldCues},'headCues',{headCues}, ...
    'angles',angles,'frameTime',t,'bearings',bearings,'distances',distances);

f = figure('Color','w','Position',[100 100 1100 760],'Visible','off');
fields = {'ILDRange_dB','ITDRange_us','PairedILDRMS_dB','PairedITDRMS_us'};
labels = {'World-fixed ILD range (dB)','World-fixed ITD range (us)', ...
    'RMS ILD difference from matched head-fixed (dB)', ...
    'RMS ITD difference from matched head-fixed (us)'};
for k = 1:4
    subplot(2,2,k);hold on;
    for d = 1:numel(distances)
        select = summary.Distance_m == distances(d);
        plot(bearings,summary.(fields{k})(select),'-o','LineWidth',1.2,'MarkerSize',4);
    end
    grid on;xlim([-150 180]);xticks(-150:60:180);
    xlabel('Bearing at turn onset (deg; positive = left)');ylabel(labels{k});
    if k == 1, legend('30 m','100 m','300 m','Location','best'); end
end
exportgraphics(f,fullfile(folder,'results','placement_cues.png'),'Resolution',150);close(f);

f = figure('Color','w','Position',[100 100 1100 650],'Visible','off');
subplot(2,2,1);hold on;
for d = 1:numel(distances)
    select = summary.Distance_m == distances(d);
    plot(bearings,summary.TranslationChange_deg(select),'-o','LineWidth',1.2);
end
yline(0,':');grid on;xlabel('Initial bearing (deg)');ylabel('Translation contribution (deg)');
legend('30 m','100 m','300 m','Location','best');
chosen = [-60, 30, 120];
for k = 1:3
    subplot(2,2,k+1);hold on;
    b = find(bearings == chosen(k));w = worldCues{2,b};h = headCues{2,b};
    plot(w.Time_s,w.ITD_us,'LineWidth',1.2);plot(h.Time_s,h.ITD_us,'--','LineWidth',1.1);
    xline(8,':');xline(12,':');xlim([7.5 12.5]);grid on;
    xlabel('Time (s)');ylabel('ITD (us)');title(sprintf('%g deg initial bearing, 100 m',chosen(k)));
    if k == 1, legend('World-fixed','Matched head-fixed','Location','best'); end
end
exportgraphics(f,fullfile(folder,'results','placement_geometry.png'),'Resolution',150);close(f);
end
