
create table res_avg as SELECT
    rho,
    epsilon,
    name,
    rot,
    AVG(gen) AS avg_gen
FROM results_rho_epsilon
GROUP BY rho, epsilon, name
ORDER BY rho, epsilon, rot;


create table res_medians_score as WITH ranked AS (
    SELECT
        rho,
        epsilon,
        name,
        rot,
        score,
        ROW_NUMBER() OVER (
            PARTITION BY rho, epsilon, name, rot
            ORDER BY score
        ) AS rn,
        COUNT(*) OVER (
            PARTITION BY rho, epsilon, name, rot
        ) AS cnt
    FROM res
)
SELECT
    rho,
    epsilon,
    name,
    rot,
    AVG(score) AS median_score
FROM ranked
WHERE rn IN ((cnt + 1) / 2, (cnt + 2) / 2)
GROUP BY rho, epsilon, name, rot
ORDER BY rho, epsilon, rot;



select rho, epsilon, sum(avg_gen) as my_sum from r_e_results GROUP by rho, epsilon ORDER by my_sum



select rho, epsilon, sum(median_gen) as my_sum from res_medians GROUP by rho, epsilon ORDER by my_sum


/*select gen median across all runs*/
SELECT gen
from results_rho_epsilon_max_ls 
ORDER by gen 
LIMIT 1 OFFSET (SELECT count(gen)/2 from (SELECT gen from results_rho_epsilon_max_ls ORDER by gen))