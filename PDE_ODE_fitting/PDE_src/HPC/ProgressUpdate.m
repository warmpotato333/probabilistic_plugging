function ProgressUpdate(total)
    persistent count;
    if isempty(count)
        count = 0;
    end
    count = count + 1;
    fprintf('Progress: %d out of %d completed.\n', count, total);
end
