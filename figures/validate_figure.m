function [ok, issues] = validate_figure(fig, spec)
% VALIDATE_FIGURE  Checks the machine-readable properties of a finished figure (size, fonts,
% panels, clipping, colors, legends, axis labels) and lists what fails.
%
% Implements : the figure checks that figure_utils.export runs on each figure it writes (Fig. 2), and that
%              fig4_tipping_input runs on Fig. 4
% Inputs     : fig    figure handle
%              spec   optional struct: nPanels (expected axes), panelLabels (e.g. {'A','B'}),
%                     yMaxExpected (largest expected y value), name (for the messages)
% Outputs    : ok       true when no check fails
%              issues   cell array of messages, one per failed check; also printed
% Seed       : none (no random draws)
%
% Checks: 1 size; 2-3 font names (Arial, Helvetica, Times, Times New Roman, Symbol) and sizes (6-12
% pt); 4 panel count; 5 panel labels; 6-8 data outside the axis limits and a y range far above
% yMaxExpected; 9 pure red or green lines; 10 legend size and empty entries; 11 axis labels; 12
% diverging colormaps not monotonic on each side; 13 axis font sizes that differ between panels;
% 14 grid lines; 15 MATLAB default line colors; 16 text that starts with "Figure" or "Fig. N".

    if nargin < 2, spec = struct(); end
    if ~isfield(spec, 'name'), spec.name = 'Figure'; end

    issues = {};

    % 1. Size: width at most 7.55 in, height at most 8.80 in, width at least 2.60 in (limits of
    % 7.5, 8.75 and 2.63 in, with a small tolerance)
    fig.Units = 'inches';
    w = fig.Position(3);
    h = fig.Position(4);

    if w > 7.55
        issues{end+1} = sprintf('[DIM] Width %.2f" exceeds the maximum 7.5"', w);
    end
    if h > 8.80
        issues{end+1} = sprintf('[DIM] Height %.2f" exceeds the maximum 8.75"', h);
    end
    if w < 2.60
        issues{end+1} = sprintf('[DIM] Width %.2f" below the minimum 2.63"', w);
    end

    % 2-3. Fonts: names among Arial, Helvetica, Times, Times New Roman and Symbol; sizes within
    % the PNAS bound of 6-12 pt after reduction, with 0.5 pt tolerance on either side
    allText = findall(fig, 'Type', 'text');
    allAxes = findall(fig, 'Type', 'axes');

    badFonts = {};
    badSizes = {};

    for i = 1:numel(allText)
        t = allText(i);
        fname = get(t, 'FontName');
        fsize = get(t, 'FontSize');
        str = get(t, 'String');

        if isempty(str), continue; end
        if iscell(str), str = strjoin(str, ' '); end

        if ~ismember(lower(fname), {'arial', 'times', 'times new roman', 'symbol', 'helvetica'})
            badFonts{end+1} = sprintf('"%s" uses %s', str, fname);
        end
        if fsize < 5.5 || fsize > 12.5
            badSizes{end+1} = sprintf('"%s" is %.0fpt', str, fsize);
        end
    end

    for i = 1:numel(allAxes)
        ax = allAxes(i);
        fname = get(ax, 'FontName');
        fsize = get(ax, 'FontSize');
        if ~ismember(lower(fname), {'arial', 'times', 'times new roman', 'symbol', 'helvetica'})
            badFonts{end+1} = sprintf('Axis uses %s', fname);
        end
        if fsize < 5.5 || fsize > 12.5
            badSizes{end+1} = sprintf('Axis tick labels are %.0fpt', fsize);
        end
    end

    if ~isempty(badFonts)
        issues{end+1} = sprintf('[FONT] Non-Arial text: %s', strjoin(unique(badFonts), '; '));
    end
    if ~isempty(badSizes)
        issues{end+1} = sprintf('[SIZE] Text outside 6-12pt: %s', strjoin(unique(badSizes), '; '));
    end

    % 4. Panel count (colorbars and legend axes left out)
    realAxes = gobjects(0);
    for i = 1:numel(allAxes)
        ax = allAxes(i);
        tag = get(ax, 'Tag');
        if contains(lower(tag), 'colorbar') || contains(lower(tag), 'legend')
            continue;
        end
        if isa(ax, 'matlab.graphics.illustration.ColorBar')
            continue;
        end
        realAxes(end+1) = ax;
    end
    nAxes = numel(realAxes);

    if isfield(spec, 'nPanels')
        if nAxes ~= spec.nPanels
            issues{end+1} = sprintf('[PANEL] Expected %d panels, found %d axes', ...
                spec.nPanels, nAxes);
        end
    end

    % 5. Panel labels: each expected label present as bold text
    if isfield(spec, 'panelLabels')
        expectedLabels = spec.panelLabels;
        foundLabels = {};
        for i = 1:numel(allText)
            t = allText(i);
            str = get(t, 'String');
            if isempty(str), continue; end
            if iscell(str), str = str{1}; end
            fw = get(t, 'FontWeight');
            if strcmp(fw, 'bold') && ismember(strtrim(str), expectedLabels)
                foundLabels{end+1} = strtrim(str);
            end
        end
        missingLabels = setdiff(expectedLabels, foundLabels);
        if ~isempty(missingLabels)
            issues{end+1} = sprintf('[LABEL] Missing panel labels: %s', ...
                strjoin(missingLabels, ', '));
        end
    end

    % 6-8. Data range: lines, scatters and vertical bars more than half an axis range outside
    % the limits, and a y maximum above 1.5 times yMaxExpected
    for i = 1:numel(realAxes)
        ax = realAxes(i);
        yl = get(ax, 'YLim');
        xl = get(ax, 'XLim');

        lines    = findall(ax, 'Type', 'line');
        bars     = findall(ax, 'Type', 'bar');
        scatters = findall(ax, 'Type', 'scatter');

        maxY = -Inf; minY = Inf;
        maxX = -Inf; minX = Inf;

        for j = 1:numel(lines)
            ydata = get(lines(j), 'YData');
            xdata = get(lines(j), 'XData');
            ydata = ydata(isfinite(ydata));
            xdata = xdata(isfinite(xdata));
            if ~isempty(ydata)
                maxY = max(maxY, max(ydata(:)));
                minY = min(minY, min(ydata(:)));
            end
            if ~isempty(xdata)
                maxX = max(maxX, max(xdata(:)));
                minX = min(minX, min(xdata(:)));
            end
        end
        for j = 1:numel(scatters)
            ydata = get(scatters(j), 'YData');
            xdata = get(scatters(j), 'XData');
            ydata = ydata(isfinite(ydata));
            xdata = xdata(isfinite(xdata));
            if ~isempty(ydata)
                maxY = max(maxY, max(ydata(:)));
                minY = min(minY, min(ydata(:)));
            end
            if ~isempty(xdata)
                maxX = max(maxX, max(xdata(:)));
                minX = min(minX, min(xdata(:)));
            end
        end
        for j = 1:numel(bars)
            barObj = bars(j);
            isHoriz = (isprop(barObj, 'Horizontal') && strcmp(get(barObj, 'Horizontal'), 'on')) || ...
                      (isprop(barObj, 'Orientation') && strcmp(get(barObj, 'Orientation'), 'horizontal'));
            if isHoriz, continue; end
            ydata = get(barObj, 'YData');
            ydata = ydata(isfinite(ydata));
            if ~isempty(ydata)
                maxY = max(maxY, max(ydata(:)));
                minY = min(minY, min(ydata(:)));
            end
        end

        % Y-data clipping (range-based, works with negative axes)
        if maxY ~= -Inf
            yRange = yl(2) - yl(1);
            if maxY > yl(2) + yRange * 0.5 || minY < yl(1) - yRange * 0.5
                issues{end+1} = sprintf('[CLIP-Y] Axis %d: Y data range [%.1f, %.1f] outside ylim [%.1f, %.1f] -- data clipped!', ...
                    i, minY, maxY, yl(1), yl(2));
            end
        end

        % X-data clipping (range-based, works with negative axes)
        if maxX ~= -Inf
            xRange = xl(2) - xl(1);
            if maxX > xl(2) + xRange * 0.5 || minX < xl(1) - xRange * 0.5
                issues{end+1} = sprintf('[CLIP-X] Axis %d: X data range [%.1f, %.1f] outside xlim [%.1f, %.1f]', ...
                    i, minX, maxX, xl(1), xl(2));
            end
        end

        % Double-scaling guard: a y value above 1.5 times spec.yMaxExpected suggests the data were scaled twice
        if isfield(spec, 'yMaxExpected') && maxY ~= -Inf
            if maxY > spec.yMaxExpected * 1.5
                issues{end+1} = sprintf('[SCALE] Axis %d: max Y data (%.1f) >> expected max (%.1f) -- possible double-scaling!', ...
                    i, maxY, spec.yMaxExpected);
            end
        end
    end

    % 9. Color-blind check: lines in pure red or pure green
    allLines = findall(fig, 'Type', 'line');
    unsafeColors = {};
    for i = 1:numel(allLines)
        c = get(allLines(i), 'Color');
        if numel(c) >= 3
            if c(1) > 0.9 && c(2) < 0.1 && c(3) < 0.1
                unsafeColors{end+1} = 'pure red [1 0 0]';
            end
            if c(1) < 0.1 && c(2) > 0.9 && c(3) < 0.1
                unsafeColors{end+1} = 'pure green [0 1 0]';
            end
        end
    end
    if ~isempty(unsafeColors)
        issues{end+1} = sprintf('[COLOR] Colorblind-unsafe colors found: %s', ...
            strjoin(unique(unsafeColors), ', '));
    end

    % 10. Legends: more than 8 entries, or an empty entry
    legends = findall(fig, 'Type', 'legend');
    for i = 1:numel(legends)
        leg = legends(i);
        nEntries = numel(leg.String);
        if nEntries > 8
            issues{end+1} = sprintf('[LEGEND] Legend has %d entries -- consider simplifying', nEntries);
        end
        for j = 1:nEntries
            if isempty(strtrim(leg.String{j}))
                issues{end+1} = sprintf('[LEGEND] Empty legend entry at position %d', j);
            end
        end
    end

    % 11. Axis labels (no x label needed with rotated ticks, no y label with horizontal bars)
    for i = 1:numel(realAxes)
        ax = realAxes(i);
        xlStr = get(get(ax, 'XLabel'), 'String');
        ylStr = get(get(ax, 'YLabel'), 'String');

        hasRotatedTicks = (get(ax, 'XTickLabelRotation') ~= 0);
        hasBarH = ~isempty(findall(ax, 'Type', 'bar', '-and', 'Horizontal', 'on'));

        if ~hasRotatedTicks && (isempty(xlStr) || (ischar(xlStr) && isempty(strtrim(xlStr))))
            issues{end+1} = sprintf('[LABEL] Axis %d missing xlabel', i);
        end
        if ~hasBarH && (isempty(ylStr) || (ischar(ylStr) && isempty(strtrim(ylStr))))
            issues{end+1} = sprintf('[LABEL] Axis %d missing ylabel', i);
        end
    end

    % 12. Colormap monotonicity (an inverted diverging colormap fails)
    for i = 1:numel(realAxes)
        ax = realAxes(i);
        imgs = findall(ax, 'Type', 'image');
        if isempty(imgs), continue; end

        cmap = colormap(ax);
        nC = size(cmap, 1);
        if nC < 16, continue; end  % skip tiny colormaps

        % Check for diverging colormap (should go dark -> light -> dark)
        % Luminance proxy: 0.299*R + 0.587*G + 0.114*B
        lum = cmap * [0.299; 0.587; 0.114];
        midIdx = round(nC / 2);
        lumStart = mean(lum(1:3));
        lumMid   = mean(lum(midIdx-1:midIdx+1));
        lumEnd   = mean(lum(end-2:end));

        % Diverging map should peak in the middle
        isDiverging = lumMid > lumStart + 0.05 && lumMid > lumEnd + 0.05;

        if isDiverging
            % Check symmetry: start and end should both be darker than middle
            % If one side is bright and middle is also bright -> not diverging properly
            firstHalf = lum(1:midIdx);
            secondHalf = lum(midIdx:end);

            % First half should be monotonically increasing (dark -> light)
            firstSmoothed = movmean(firstHalf, 5);
            firstMonotonic = all(diff(firstSmoothed) >= -0.03);

            % Second half should be monotonically decreasing (light -> dark)
            secondSmoothed = movmean(secondHalf, 5);
            secondMonotonic = all(diff(secondSmoothed) <= 0.03);

            if ~firstMonotonic || ~secondMonotonic
                issues{end+1} = sprintf('[CMAP] Axis %d: diverging colormap not monotonic (possible inversion)', i);
            end
        end
    end

    % 13. Axis font sizes differing between panels by more than 0.5 pt
    if numel(realAxes) > 1
        axFontSizes = zeros(numel(realAxes), 1);
        for i = 1:numel(realAxes)
            axFontSizes(i) = get(realAxes(i), 'FontSize');
        end
        if range(axFontSizes) > 0.5
            issues{end+1} = sprintf('[CONSIST] Axes font sizes vary: %s pt', ...
                strjoin(arrayfun(@(s) sprintf('%.0f', s), axFontSizes, 'UniformOutput', false), ', '));
        end
    end

    % 14. Grid lines (Rule 8 of Rougier et al. 2014, "Ten simple rules for better figures")
    for i = 1:numel(realAxes)
        ax = realAxes(i);
        if strcmp(get(ax, 'XGrid'), 'on') || strcmp(get(ax, 'YGrid'), 'on')
            issues{end+1} = sprintf('[GRID] Axis %d has grid lines enabled (consider removing for cleaner look)', i);
        end
    end

    % 15. MATLAB default line colors (the first four of its default order), which are not in
    % the Okabe-Ito palette
    matlabDefaults = [
        0, 0.4470, 0.7410;   % default blue
        0.8500, 0.3250, 0.0980;  % default orange
        0.9290, 0.6940, 0.1250;  % default yellow
        0.4940, 0.1840, 0.5560;  % default purple
    ];
    defaultFound = {};
    allPlotLines = findall(fig, 'Type', 'line');
    for i = 1:numel(allPlotLines)
        c = get(allPlotLines(i), 'Color');
        if numel(c) < 3, continue; end
        for d = 1:size(matlabDefaults, 1)
            if max(abs(c(1:3) - matlabDefaults(d,:))) < 0.01
                defaultFound{end+1} = sprintf('[%.3f %.3f %.3f]', matlabDefaults(d,:));
                break;
            end
        end
    end
    if ~isempty(defaultFound)
        issues{end+1} = sprintf('[DEFAULT] MATLAB default colors detected: %s (use Okabe-Ito palette)', ...
            strjoin(unique(defaultFound), ', '));
    end

    % 16. Embedded figure titles: figure numbers, titles and captions belong in the caption
    for i = 1:numel(allText)
        t = allText(i);
        str = get(t, 'String');
        if isempty(str), continue; end
        if iscell(str), str = strjoin(str, ' '); end
        strLower = lower(strtrim(str));
        % Check for embedded "Figure N" or "Fig N" text (not panel labels)
        if ~isempty(regexp(strLower, '^(figure|fig\.?\s*\d)', 'once'))
            issues{end+1} = sprintf('[EMBED] Embedded figure title detected: "%s" (move to caption)', str);
        end
    end

    % Result: ok when no check failed; the messages are printed
    nErrors = numel(issues);
    ok = (nErrors == 0);

    if ok
        fprintf('  VALIDATE %s: PASS (0 issues)\n', spec.name);
    else
        fprintf('  VALIDATE %s: %d issue(s):\n', spec.name, nErrors);
        for i = 1:nErrors
            fprintf('    %s\n', issues{i});
        end
    end
end
