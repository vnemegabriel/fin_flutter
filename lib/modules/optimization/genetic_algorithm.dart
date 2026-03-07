import 'dart:math';
import 'nelder_mead.dart';

/// Genetic Algorithm with real-coded chromosomes.
///
/// Uses Simulated Binary Crossover (SBX) and polynomial mutation.
/// Suitable for multi-modal problems and larger variable spaces.
class GeneticAlgorithm {
  final int populationSize;
  final int maxGenerations;
  final double crossoverRate;
  final double mutationRate;
  final double distributionIndexCrossover; // eta_c for SBX (default 20)
  final double distributionIndexMutation;  // eta_m for polynomial mutation (default 20)
  final int eliteCount;

  final Random _rng;

  GeneticAlgorithm({
    this.populationSize = 50,
    this.maxGenerations = 200,
    this.crossoverRate = 0.9,
    this.mutationRate = 0.05,
    this.distributionIndexCrossover = 20.0,
    this.distributionIndexMutation = 20.0,
    this.eliteCount = 2,
    int seed = 42,
  }) : _rng = Random(seed);

  /// Optimize objective function.
  ///
  /// Returns stream of iteration results.
  Iterable<OptimizationIteration> optimize({
    required ObjectiveFunction objective,
    required int nVariables,
  }) sync* {
    // Initialize population
    var population = List.generate(
        populationSize, (_) => List.generate(nVariables, (_) => _rng.nextDouble()));
    var fitness = population.map(objective).toList();

    for (var gen = 0; gen < maxGenerations; gen++) {
      // Sort by fitness
      final indices = List.generate(populationSize, (i) => i)
        ..sort((a, b) => fitness[a].compareTo(fitness[b]));

      yield OptimizationIteration(
        iteration: gen,
        bestValue: fitness[indices[0]],
        bestPoint: population[indices[0]],
        converged: false,
      );

      // New population
      final newPop = <List<double>>[];

      // Elitism: preserve top individuals
      for (var i = 0; i < eliteCount; i++) {
        newPop.add(List<double>.from(population[indices[i]]));
      }

      // Fill rest with crossover + mutation
      while (newPop.length < populationSize) {
        final p1 = _tournamentSelect(population, fitness);
        final p2 = _tournamentSelect(population, fitness);

        List<double> child1, child2;
        if (_rng.nextDouble() < crossoverRate) {
          (child1, child2) = _sbxCrossover(p1, p2);
        } else {
          child1 = List<double>.from(p1);
          child2 = List<double>.from(p2);
        }

        _polynomialMutation(child1);
        _polynomialMutation(child2);

        newPop.add(child1);
        if (newPop.length < populationSize) newPop.add(child2);
      }

      population = newPop;
      fitness = population.map(objective).toList();
    }

    // Final result
    final bestIdx = List.generate(populationSize, (i) => i)
        .reduce((a, b) => fitness[a] < fitness[b] ? a : b);
    yield OptimizationIteration(
      iteration: maxGenerations,
      bestValue: fitness[bestIdx],
      bestPoint: population[bestIdx],
      converged: true,
    );
  }

  /// Tournament selection (size 3).
  List<double> _tournamentSelect(
      List<List<double>> population, List<double> fitness) {
    var best = _rng.nextInt(populationSize);
    for (var i = 0; i < 2; i++) {
      final candidate = _rng.nextInt(populationSize);
      if (fitness[candidate] < fitness[best]) best = candidate;
    }
    return population[best];
  }

  /// Simulated Binary Crossover (SBX).
  (List<double>, List<double>) _sbxCrossover(
      List<double> p1, List<double> p2) {
    final n = p1.length;
    final c1 = List<double>.from(p1);
    final c2 = List<double>.from(p2);

    for (var i = 0; i < n; i++) {
      if (_rng.nextDouble() > 0.5) continue;
      final u = _rng.nextDouble();
      double beta;
      if (u <= 0.5) {
        beta = pow(2.0 * u, 1.0 / (distributionIndexCrossover + 1.0)).toDouble();
      } else {
        beta = pow(1.0 / (2.0 * (1.0 - u)), 1.0 / (distributionIndexCrossover + 1.0))
            .toDouble();
      }
      c1[i] = 0.5 * ((1 + beta) * p1[i] + (1 - beta) * p2[i]);
      c2[i] = 0.5 * ((1 - beta) * p1[i] + (1 + beta) * p2[i]);
      c1[i] = c1[i].clamp(0.0, 1.0);
      c2[i] = c2[i].clamp(0.0, 1.0);
    }
    return (c1, c2);
  }

  /// Polynomial mutation.
  void _polynomialMutation(List<double> x) {
    for (var i = 0; i < x.length; i++) {
      if (_rng.nextDouble() > mutationRate) continue;
      final u = _rng.nextDouble();
      double delta;
      if (u < 0.5) {
        delta = pow(2.0 * u, 1.0 / (distributionIndexMutation + 1.0)).toDouble() - 1.0;
      } else {
        delta = 1.0 - pow(2.0 * (1.0 - u), 1.0 / (distributionIndexMutation + 1.0))
            .toDouble();
      }
      x[i] = (x[i] + delta).clamp(0.0, 1.0);
    }
  }
}
