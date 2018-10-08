#include <iostream>    // std::cout
#include <string>      // std::string
#include <memory>      // std::unique_ptr
#include <algorithm>   // std::generate_n
#include <numeric>     // std::inner_product
#include <chrono>      // std::chrono
#include <iomanip>     // std::setw
#include <random>      // std::uniform_distribution
#include <type_traits> // std::result_of, std::invoke_result
#ifdef _OPENMP
#include <omp.h>       // omp_get_max_threads
#include <merge_sort.h>
#endif

#include <mass.h>

static constexpr int rng_seed    = 123;
static constexpr float inf       = std::numeric_limits<float>::infinity();
static constexpr float fit_limit = 0.f;


auto printProgress = [](const float &now, const int &total, const std::chrono::high_resolution_clock::time_point &start_time)
                      {
                        float perc = now / total;
                        int lpad = static_cast<int>(perc * PBWIDTH);
                        std::cout << "\rOptimization progress:"
                                  << std::right << std::setw(5) << std::setprecision(3) << perc * 100.f << "%  ["
                                  << std::left  << std::setw(PBWIDTH - 1) << std::string(lpad, '|') << "] "
                                  << std::right << std::setw(5) << int(now) << "/" << std::setprecision(5) << total
                                  << " [" << std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - start_time).count() << " sec]";
                        std::cout << std::flush;
                      };


template<typename coords, typename laplacian, typename Func>
auto par_gen(const coords &target,
             const laplacian &lap,
             const int &Nnodes,
             const int &n_population,
             const int &max_iter,
             Func fit,
             float elit_rate = .1f,
             float mutation_rate = .3f,
             std::size_t seed = 0,
             int verbose = 1, // logfile
             int nth = 4)
{
  using res_t = typename std::result_of<Func(const coords &, const coords &)>::type; // since c++17
#ifdef _OPENMP
  nth -= nth % 2;
  const int diff_size = n_population % nth,
            size      = diff_size ? n_population - diff_size : n_population;
#endif
  const int elite = n_population * elit_rate;
  int iteration = 0,
      best = 0;
  std::unique_ptr<masses[]> population(new masses[n_population]),
                            new_gen(   new masses[n_population]);
  std::unique_ptr<int[]>    rank(      new int[n_population]);
  std::unique_ptr<res_t[]> fitness (   new res_t[n_population]);

  auto eigen_lap = Eigen::Map<Eigen::Matrix<float,
                              Eigen::Dynamic,
                              Eigen::Dynamic,
                              Eigen::RowMajor>>
                              (&lap[0], N, N);

  std::mt19937 engine(seed);
  std::uniform_real_distribution<float> rng;
  auto start_time = std::chrono::high_resolution_clock::now();


#ifdef _OPENMP
#pragma omp parallel num_threads(nth)
  {
#endif

#ifdef _OPENMP
#pragma omp for
  for(int i = 0; i < n_population; ++i) population[i] = masses(Nnodes, rng(engine));
#else
  std::generate_n(population.get(), n_population, [&](){return masses(Nnodes, rng(engine));});
#endif

  while(iteration < max_iter)
  {
#ifdef _OPENMP
#pragma omp for
    for(int i = 0; i < n_population; ++i)
    {
      // argsort variables
      fitness[i] = fit(eigen_lap, population[i], target);
      rank[i]    = i;
    }
#else
    std::transform(population.get(), population.get() + n_population,
                   fitness.get(),
                   [&](const genome &pop)
                   {
                    return fit(eigen_lap, pop, target);
                  });
    std::iota(rank.get(), rank.get() + n_population, 0);
#endif

#ifdef _OPENMP
#pragma omp single
    {
      mergeargsort_parallel_omp(rank.get(), fitness.get(), 0, size, nth, [&](const int &a1, const int &a2){return fitness[a1] < fitness[a2];});
      if (diff_size)
      {
        std::sort(rank.get() + size, rank.get() + n_population, [&](const int &a1, const int &a2){return fitness[a1] < fitness[a2];});
        std::inplace_merge(rank.get(), rank.get() + size, rank.get() + n_population, [&](const int &a1, const int &a2){return fitness[a1] < fitness[a2];});
      }
    }
#else
    std::sort(rank.get(), rank.get() + n_population, [&](const int &a1, const int &a2){return fitness[a1] < fitness[a2];});
#endif

    best = rank[0];

#ifdef _OPENMP
#pragma omp single nowait
    {
#endif
      switch(verbose)
      {
        case 1: printProgress(match(population[best], target), Nnodes, start_time);
        break;
        case 2:
        {
          std::cout << "iter: "
                    << std::setw(5)          << iteration        << "   "
          //          << std::setw(LENGTH + 3) << population[best] << "   "
                    << std::setw(5)          << fitness[best]    << "   "
          //          << std::setw(LENGTH + 3) << target           << "   "
          //          << std::setw(3)          << match(population[best], target)
                    << std::endl;
        } break;
        default: break;
      }
#ifdef _OPENMP
    } // close single section
#endif
    if (fitness[best] == 0) break;

#ifdef _OPENMP
#pragma omp for
#endif
    for(int i = 0; i < n_population; ++i)
    {
      // crossover
      new_gen[i] = ( rank[i] < elite ) ? population[rank[i]]
                                       : population[rank[int(rng(engine) * elite)]] + population[rank[int(rng(engine) * elite)]];
      // random mutation
      if (rng(engine) < mutation_rate) !new_gen[i];
    }

#ifdef _OPENMP
#pragma omp for
    for(int i = 0; i < n_population; ++i) population[i] = new_gen[i];
#else
    std::copy_n(new_gen.get(), n_population, population.get());
#endif

#ifdef _OPENMP
#pragma omp single
#endif
    ++iteration;

  } // end while

#ifdef _OPENMP
  } // end parallel section
#endif

  if (verbose == 1) std::cout << std::endl;

  return population[best];
}


