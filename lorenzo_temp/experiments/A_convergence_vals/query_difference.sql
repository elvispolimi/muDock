SELECT SUM(difference) AS total_difference
FROM (
    SELECT
        autostop.score - noautostop.score AS difference
    FROM autostop
    JOIN noautostop
        ON autostop.name = noautostop.name
        AND autostop.seed = noautostop.seed
);