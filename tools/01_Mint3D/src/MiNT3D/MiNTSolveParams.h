#pragma once

/*
 * This struct is used to store all of the parameters needed to specify a solve at a specific outer iteration.
 */

#include <iomanip>   // Required for std::setprecision and std::scientific
#include <iostream>

// #include "../mint3D_hook/Serialization.h"

namespace MiNT3D {

enum class BoundaryCondition { Free = 0, SDF, Poisson };

enum class BoundaryHardConstraints { NoHardConstraints = 0, Normal, Frame };

enum class InitState { Exact = 0, Random };

enum class ConnectionForCombing { AsIntegrableAsPossible = 0, SmoothnessUsesAsSmoothAsPossible };

enum class FitType { Moments = 0, Direction, Kruskal };

class SolveParams {
   public:
    // I want to add a section which initializes SolveParams objects with different choices of default values

    // constructor with default values
    SolveParams() {
        // use call to resetDefault() to set default values
        *this = resetDefault();
    }

    SolveParams(const std::string& filepath) {
        *this = resetDefault();
        // load parameters from file
        // Serialization::deserializeSolveParams(filepath, *this);
    }

    void setFitSolve() {
        // Additional energy weights from GUI
        this->w_smooth_combed = 1.;          // 10
        this->w_smooth_combed_precon = 0.;   // 1e5
        this->w_smooth_sym = 1.;
        this->w_smooth_sym_precon = 0;   //  100;
        this->w_int_combed = 1.;
        this->w_int_combed_fixed_weight = 1.;   // 1e-3;
        this->w_int_sym = 1.;
        this->w_int_sym_fixed_weight = 0.;
        this->w_orthog = 5e1;
        this->w_orthog_precon = 0.;
        this->w_orthog_constraint_weight = 0.;
        this->w_unit = 1e-6;   // also ok 8e-3
        this->w_unit_precon = 0.;
        this->w_fit = 1e6;   // also ok 1e-3
        this->w_fit_precon = 0;
        this->w_self_align = 0;
        this->w_viscocity = 1e-4;
        this->w_feat_align = 0;
        this->w_init_scale = 1.;   //  1e-4

        this->w_max_lambda = 1e17;

        // Booleans from GUI
        this->b_smooth_combed = false;
        this->b_smooth_sym = true;
        this->b_smooth_asap_combing = false;
        this->b_int_combed = false;
        this->b_int_sym = true;
        this->b_orthog = false;
        this->b_unit = true;
        this->b_fit = true;
        this->b_self_align = false;
        this->b_viscosity = true;
        this->b_feat_align = false;
        this->b_init_as_odeco = false;

        this->fit_type = FitType::Direction;

        this->boundary_hard_constraints =
            BoundaryHardConstraints::NoHardConstraints;   //  BoundaryHardConstraints::Normal;
    }

    void setAlignWithNormalBoundarySolve() {
        // Additional energy weights from GUI
        this->w_smooth_combed = 1.;          // 10
        this->w_smooth_combed_precon = 0.;   // 1e5
        this->w_smooth_sym = 1.;
        this->w_smooth_sym_precon = 0;   //  100;
        this->w_int_combed = 1.;
        this->w_int_combed_fixed_weight = 1.;   // 1e-3;
        this->w_int_sym = 1.;
        this->w_int_sym_fixed_weight = 0.;
        this->w_orthog = 5e1;
        this->w_orthog_precon = 0.;
        this->w_orthog_constraint_weight = 0.;
        this->w_unit = 1e-4;   // also ok 8e-3
        this->w_unit_precon = 0.;
        this->w_fit = 5e-2;   // also ok 1e-3
        this->w_fit_precon = 0;
        this->w_self_align = 0;
        this->w_viscocity = 0;
        this->w_feat_align = 0;
        this->w_init_scale = 1.;   //  1e-4

        this->w_max_lambda = 1e17;

        // Booleans from GUI
        this->b_smooth_combed = false;
        this->b_smooth_sym = true;
        this->b_smooth_asap_combing = false;
        this->b_int_combed = false;
        this->b_int_sym = true;
        this->b_orthog = false;
        this->b_unit = true;
        this->b_fit = true;
        this->b_self_align = false;
        this->b_viscosity = false;
        this->b_feat_align = false;
        this->b_init_as_odeco = false;

        this->fit_type = FitType::Direction;

        this->boundary_hard_constraints = BoundaryHardConstraints::Normal;   //  BoundaryHardConstraints::Normal;
    }

    void setSoftAlignSolve() {
        // Additional energy weights from GUI
        this->w_smooth_combed = 1.;          // 10
        this->w_smooth_combed_precon = 0.;   // 1e5
        this->w_smooth_sym = 1.;
        this->w_smooth_sym_precon = 0;   //  100;
        this->w_int_combed = 1.;
        this->w_int_combed_fixed_weight = 1.;   // 1e-3;
        this->w_int_sym = 1e-5;
        this->w_int_sym_fixed_weight = 0.;
        this->w_orthog = 100;
        this->w_orthog_precon = 0.;
        this->w_orthog_constraint_weight = 0.;
        this->w_unit = 2.5;   //  5e-3;   // also ok 8e-3
        this->w_unit_precon = 0.;
        this->w_fit = .01;   // 5e-4;   // also ok 1e-3
        this->w_fit_precon = 0;
        this->w_self_align = 0;
        this->w_viscocity = this->w_int_sym * 1.e-3;
        this->w_feat_align = 0;

        this->w_init_scale = 1.;   //  1e-4

        this->w_max_lambda = 1e17;

        // Booleans from GUI
        this->b_smooth_combed = false;
        this->b_smooth_sym = true;
        this->b_smooth_asap_combing = false;
        this->b_int_combed = false;
        this->b_int_sym = true;
        this->b_orthog = true;
        this->b_unit = true;
        this->b_fit = true;
        this->b_self_align = false;
        this->b_viscosity = true;
        this->b_feat_align = false;
        this->b_init_as_odeco = false;

        this->fit_type = FitType::Direction;

        this->boundary_hard_constraints =
            BoundaryHardConstraints::NoHardConstraints;   //  BoundaryHardConstraints::Normal;
    }

    // the weights here don't matter, they get overwritten
    void setOdecoSolve() {
        // Additional energy weights from GUI

        this->w_smooth_combed = 1.;          // 10
        this->w_smooth_combed_precon = 0.;   // 1e5
        this->w_smooth_sym = 1e-2;        // 1e0;            // 1e-12;  / GL config           // 1e-4;           // .1;
        this->w_smooth_sym_precon = 0.;   // 1e11;    // 1e3;
        this->w_int_combed = 1.;
        this->w_int_combed_fixed_weight = 0.;   // 1e-3;
        this->w_int_sym = 1.;
        this->w_int_sym_fixed_weight = 0.;
        this->w_orthog = 50.;        // 1e-4;
        this->w_orthog_precon = 0;   // 1e12;
        this->w_orthog_constraint_weight = 0.;
        this->w_unit = 1e0;         // 1e-8;        //.0005;       // 1e3;   // 1e3;
        this->w_unit_precon = 0.;   // 1e6;   // 1e3;   // 1e8;
        this->w_fit = 0;
        this->w_fit_precon = 0;
        this->w_self_align = 0;
        this->w_viscocity = 0;
        this->w_feat_align = 0;

        this->w_init_scale = .1;   //  1e-4

        // Booleans from GUI
        this->b_smooth_combed = false;
        this->b_smooth_sym = true;
        this->b_smooth_asap_combing = false;
        this->b_int_combed = false;
        this->b_int_sym = false;
        this->b_orthog = true;   // true;
        this->b_unit = true;     // true;
        this->b_fit = false;     // this is hacky, seperate these two things
        this->b_self_align = false;
        this->b_viscosity = false;
        this->b_feat_align = false;
        this->b_init_as_odeco = false;

        this->w_max_lambda = 8;   // 10. * this->w_orthog / this->w_int_sym;
        // 1e10;   //(this->w_orthog / 1.e0) * 1. / this->w_int_sym;

        this->boundary_hard_constraints = BoundaryHardConstraints::Normal;
    }

    void setSoftOrthogSolve() {
        // Additional energy weights from GUI

        this->w_outer_step = 1.414;   // std::pow(2, .25);   // 414;   // 1.414;   // 1.1;
        this->inner_iter_max_steps = 50;

        this->w_smooth_combed = 1.;          // 10
        this->w_smooth_combed_precon = 0.;   // 1e5
        this->w_smooth_sym = 1.;             // 1e-12;  / GL config           // 1e-4;           // .1;
        this->w_smooth_sym_precon = 0;       // 1e10;    // 1e11;    // 1e3;
        this->w_int_combed = 0.;
        this->w_int_combed_fixed_weight = 0.;   // 1e-3;
        this->w_int_sym = 1.e-8;
        this->w_int_sym_fixed_weight = 0.;
        this->w_orthog = .1;
        this->w_orthog_precon = 0;   // 1e9;   // 1e12;
        this->w_orthog_constraint_weight = 0.;
        // can change to get diff results 1e3-1e-2 is a good range to explore
        this->w_unit = 1.e-5;       // 1e3;   // 1e3;
        this->w_unit_precon = 0.;   // 1e6;   // 1e3;   // 1e8;
        this->w_fit = 0;
        this->w_fit_precon = 0;
        this->w_self_align =
            this->w_int_sym *
            1.e-7;   // 1. / this->w_orthog * 1e2;   // this will be 1e2 when int penalty equals orthog penalty
        this->w_viscocity = 1.e-9;
        this->w_feat_align = 0.;
        this->w_init_scale = 1e-6;   // 5.;   //  1e-4

        // Booleans from GUI
        this->b_smooth_combed = false;
        this->b_smooth_sym = true;
        this->b_smooth_asap_combing = false;
        this->b_int_combed = false;
        this->b_int_sym = true;
        this->b_orthog = true;
        this->b_unit = true;
        this->b_fit = false;   // this is hacky, seperate these two things
        this->b_self_align = false;
        this->b_viscosity = true;
        this->b_feat_align = true;
        this->b_init_as_odeco = false;

        // this->b_use_kruskal_tensors_for_sym_smoothness = true;

        this->w_max_lambda = 1e17;   // 10. * this->w_orthog / this->w_int_sym;
        // 1e10;   //(this->w_orthog / 1.e0) * 1. / this->w_int_sym;

        this->boundary_hard_constraints = BoundaryHardConstraints::Normal;

        /*
        this->w_smooth_combed = .1;          // 10
        this->w_smooth_combed_precon = 0.;   // 1e5
        this->w_smooth_sym = 1.;             // 1e-4;           // .1;
        this->w_smooth_sym_precon = 1e16;    // 1e3;
        this->w_int_combed = 0.;
        this->w_int_combed_fixed_weight = 1e-4;   // 1e-3;
        this->w_int_sym = 1e-2;
        this->w_int_sym_fixed_weight = 0.;
        this->w_orthog = 1e14;
        this->w_orthog_precon = 0.;
        this->w_orthog_constraint_weight = 0.;
        this->w_unit = 1e1;   /// not actually soft but w/e this is stable
        this->w_unit_precon = 0.;
        this->w_fit = 0;
        this->w_fit_precon = 0;
        this->w_init_scale = 1.;   //  1e-4

        // Booleans from GUI
        this->b_smooth_combed = false;
        this->b_smooth_sym = true;
        this->b_smooth_asap_combing = false;
        this->b_int_combed = false;
        this->b_int_sym = true;
        this->b_orthog = true;
        this->b_unit = true;
        this->b_fit = false;

        this->w_max_lambda = 1e17;   //(this->w_orthog / 1.e0) * 1. / this->w_int_sym;



*/

        // // Additional energy weights from GUI
        // this->w_smooth_combed = 1e-3;        // 10
        // this->w_smooth_combed_precon = 0.;   // 1e5
        // this->w_smooth_sym = 1e0;
        // this->w_smooth_sym_precon = 1e5;
        // this->w_int_combed = 0.;
        // this->w_int_combed_fixed_weight = 1e-3;
        // this->w_int_sym = 1e-2;
        // this->w_int_sym_fixed_weight = 0.;
        // this->w_orthog = 1e1;
        // this->w_orthog_precon = 1e4;
        // this->w_unit = 5e-3;
        // this->w_unit_precon = 1e3;
        // this->w_init_scale = 1.;   //  1e-4

        // this->w_max_lambda = 1e10;

        // // Booleans from GUI
        // this->b_smooth_combed = true;
        // this->b_smooth_sym = true;
        // this->b_smooth_asap_combing = false;
        // this->b_int_combed = true;
        // this->b_int_sym = true;
        // this->b_orthog = true;
        // this->b_unit = true;
    }

    SolveParams resetDefault() {
        SolveParams new_params = *this;

        // Solver settings

        new_params.inner_iter_max_steps = 50;   // 50
        new_params.grad_tol = 1e-8;
        new_params.xTol = 0;       // do not consider change in x as a stopping criterion
        new_params.fTol = 1e-14;   // use very small change in objective as one stopping criterion
        new_params.reg = 1e-8;     // initial magnitude of the diagonal regularization term in the newton solver

        new_params.total_step = 0;
        new_params.outer_step = 0;

        // This controls what tensor representation to use for the integrability term (krushkal vs moment)
        // we use krushkal for smoothness and moment for integrability
        // In theory it is ok to use kruskal for moment as long as rank doesn't degenerate.

        new_params.b_use_kruskal_tensors_for_sym_smoothness = false;

        new_params.w_outer_step = 1.414;   // 1.1;
        new_params.lambda_penalty = 1.;

        // Additional energy weights from GUI
        new_params.w_smooth_combed = 1e-3;        // 10
        new_params.w_smooth_combed_precon = 0.;   // 1e5
        new_params.w_smooth_sym = 1e0;
        new_params.w_smooth_sym_precon = 1e5;
        new_params.w_int_combed = 0.;
        new_params.w_int_combed_fixed_weight = 1e-3;
        new_params.w_int_sym = 1e-2;
        new_params.w_int_sym_fixed_weight = 0.;
        new_params.w_orthog = 1e1;
        new_params.w_orthog_precon = 1e4;
        new_params.w_orthog_constraint_weight = 0.;
        new_params.w_unit = 1e2;
        new_params.w_unit_precon = 1e3;
        new_params.w_fit = 0;
        new_params.w_fit_precon = 0;
        new_params.w_self_align = 0;
        new_params.w_viscocity = 1e-5;
        new_params.w_feat_align = 1e4;
        new_params.w_init_scale = 1.;   //  1e-4

        // Booleans from GUI
        new_params.b_smooth_combed = false;
        new_params.b_smooth_sym = true;
        new_params.b_smooth_asap_combing = false;
        new_params.b_int_combed = false;
        new_params.b_int_sym = true;
        new_params.b_orthog = false;
        new_params.b_unit = false;
        new_params.b_fit = false;
        new_params.b_self_align = false;
        new_params.b_viscosity = true;
        new_params.b_feat_align = true;
        new_params.b_init_as_odeco = false;

        // Boundary conditions and selected options
        new_params.boundary_condition = BoundaryCondition::SDF;   // BoundaryCondition::Poisson;
        new_params.boundary_hard_constraints = BoundaryHardConstraints::Normal;
        new_params.init_state = InitState::Exact;   // InitState::Random;
        new_params.connection_for_combing = ConnectionForCombing::AsIntegrableAsPossible;

        new_params.w_bound = 0.;
        new_params.w_bound_rescale = 1.;
        new_params.w_mint = 1.;   // 1e-12;   // 1e-8;
        new_params.w_max_lambda = 1e17;
        new_params.w_smooth = 1;
        new_params.w_scaled_jacobian = 0.;
        new_params.w_unit_barrier = 1. * 0.;
        new_params.w_unit_norm = 1.;
        // new_params.w_fit = 0.;

        new_params.smoothness_energy = 0;
        new_params.unit_energy = 0;
        new_params.mint_energy = 0;
        new_params.primal_integrability_energy = 0;
        new_params.asap_combed_smoothness_energy = 0;
        new_params.combed_smoothness_energy = 0;
        new_params.scaled_jacobian_energy = 0;
        new_params.unit_barrier_energy = 0;
        new_params.total_energy = 0;
        new_params.fit_energy = 0;

        return new_params;
    }

    // I want a constructor that initializes the SolveParams object with a specific set of values from file

    SolveParams setDefaultWithScale(double scale) {
        SolveParams new_params = resetDefault();
        // new_params.w_global_rescale_param = scale;
        new_params.setGlobalScale(scale);
        return new_params;
    }

    // Solver settings

    int inner_iter_max_steps = 50;
    double grad_tol = 1e-8;
    double xTol = 0;
    // double fTol = 1e-8;
    // double fTol = 1e-10;

    double fTol = 1e-14;   // use this is mint term uses unscaled moment, not kruskals
    double reg = 1e-8;

    int total_step = 0;
    int outer_step = 0;

    bool b_use_kruskal_tensors_for_sym_smoothness = false;
    bool b_headless_mode = false;

    double w_outer_step = 1.414;
    double lambda_penalty = 1.;

    // Energy Lambda weights

    double w_bound = 1.;   // uniform scale applied to all boundary vecs.

    // relative scale of frame vectors.
    // On boundary can use this to control anisotropy of normal direction
    // vs regularity of vectors in the tangent space.
    // double w_bound_vec0 = 1. / 3.;
    // double w_bound_vec1 = 1. / 3.;
    // double w_bound_vec2 = 1. / 3.;

    double w_bound_rescale = 1.;

    // Current weight of moment integrability term
    double w_mint = 1e-4;
    // double w_mint = 1e15;

    // Threshold for when to terminate the outer loop
    double w_max_lambda = 1e17;   // 1e17;

    // Additional energy weights from GUI
    double w_smooth_combed = 1.0;
    double w_smooth_combed_precon = 1.0;
    double w_smooth_sym = 1.0;
    double w_smooth_sym_precon = 1.0;
    double w_int_combed = 1.0;
    double w_int_combed_fixed_weight = 1.0;
    double w_int_sym = 1.0;
    double w_int_sym_fixed_weight = 1.0;
    double w_orthog = 1.0;
    double w_orthog_precon = 1.0;
    double w_orthog_constraint_weight = 0.;
    double w_unit = 1.0;
    double w_unit_precon = 1.0;
    double w_fit = 0;
    double w_fit_precon = 0;
    double w_self_align = 0;
    // double w_self_align_const_reg = 0;
    double w_viscocity = 0;
    double w_feat_align = 0;
    double w_init_scale = 1e-6;

    // Booleans from GUI
    bool b_smooth_combed = false;
    bool b_smooth_sym = false;
    bool b_smooth_asap_combing = false;
    bool b_int_combed = false;
    bool b_int_sym = false;
    bool b_orthog = false;
    bool b_unit = false;
    bool b_fit = false;
    bool b_self_align = false;
    bool b_viscosity = false;
    bool b_feat_align = false;
    bool b_init_as_odeco = false;

    // Boundary conditions and selected options
    BoundaryCondition boundary_condition = BoundaryCondition::Poisson;
    BoundaryHardConstraints boundary_hard_constraints = BoundaryHardConstraints::Normal;
    InitState init_state = InitState::Exact;
    ConnectionForCombing connection_for_combing = ConnectionForCombing::AsIntegrableAsPossible;
    FitType fit_type = FitType::Direction;

    ////////////// Remove these!

    // set to 1 in "natural" units.
    // drive to zero slowly to recover the "ginzburg-landau" limit.
    double w_smooth = 1;

    // Threshold for when to terminate the outer loop
    double w_scaled_jacobian = 100;   // 100;

    // Threshold for when to terminate the outer loop
    double w_unit_barrier = 0.1;

    // The two parts of the orthonormality constraint
    double w_unit_norm = 1.;

    // This is the penalty weight of the fitting term.
    // should only be turned on if this is more than zero.
    // double w_fit = 0.;   //

    // state from step
    double smoothness_energy = 0;
    double unit_energy = 0;
    double mint_energy = 0;
    double primal_integrability_energy = 0;
    double asap_combed_smoothness_energy = 0;
    double aiap_combed_smoothness_energy = 0;
    double combed_smoothness_energy = 0;
    double scaled_jacobian_energy = 0;
    double unit_barrier_energy = 0;
    double total_energy = 0;
    double fit_energy = 0;

    SolveParams setGlobalScale(double new_targ_scale) {
        SolveParams new_weights = *this;
        double prev_global_scale = new_weights.w_global_rescale_param;
        new_weights.w_global_rescale_param = new_targ_scale;

        double scale = new_targ_scale / prev_global_scale;

        // new_weights.w_smooth_combed = w_smooth_combed * scale;
        // new_weights.w_smooth_combed_precon = w_smooth_combed_precon * scale;
        // new_weights.w_smooth_sym = w_smooth_sym * scale;
        // new_weights.w_smooth_sym_precon = w_smooth_sym_precon * scale;
        // new_weights.w_int_combed = w_int_combed * scale;
        // new_weights.w_int_combed_precon = w_int_combed_precon * scale;
        // new_weights.w_int_sym = w_int_sym * scale;
        // new_weights.w_int_sym_precon = w_int_sym_precon * scale;
        // new_weights.w_orthog = w_orthog * scale;
        // new_weights.w_orthog_precon = w_orthog_precon * scale;
        // new_weights.w_unit = w_unit * scale;
        // new_weights.w_unit_precon = w_unit_precon * scale;

        new_weights.w_bound = w_bound * scale;
        new_weights.w_mint = w_mint * scale;
        new_weights.w_smooth = w_smooth * scale;
        new_weights.w_scaled_jacobian = w_scaled_jacobian * scale;
        new_weights.w_unit_barrier = w_unit_barrier * scale;
        new_weights.w_unit_norm = w_unit_norm * scale;
        // new_weights.w_fit = w_fit * scale;

        // std::cout << "Global scale: " << new_targ_scale << " ";

        return new_weights;
    }

    std::string to_string() {
        auto format = [](const std::string& name, auto value) {
            std::ostringstream out;
            out << std::scientific << std::setprecision(5) << value;
            return name + ": " + out.str();
        };

        std::string ret = "no SCALE TERM \n\n";

        double scale = 1. / w_global_rescale_param;

        ret += format("w_bound", w_bound * scale);
        ret += ", " + format("w_mint", w_mint * scale);
        ret += ", " + format("w_smooth", w_smooth * scale);
        ret += ", " + format("w_scaled_jacobian", w_scaled_jacobian * scale);
        ret += ", " + format("w_unit_norm", w_unit_norm * scale) + "\n";
        ret += ", " + format("w_unit_barrier", w_unit_barrier * scale);
        ret += ", " + format("w_global_rescale_param", w_global_rescale_param * scale);
        // ret += ", " + format("w_fit", w_fit * scale);

        // ret += "\n\nwith SCALE TERM: ";
        // ret += format("w_global_rescale_param", w_global_rescale_param) + "\n\n";

        // ret += format("w_bound", w_bound);
        // ret += ", " + format("w_bound_viscosity", w_bound_viscosity);
        // ret += ", " + format("w_mint", w_mint);
        // ret += ", " + format("w_smooth", w_smooth);
        // ret += ", " + format("w_scaled_jacobian", w_scaled_jacobian);
        // ret += ", " + format("w_unit_norm", w_unit_norm) + "\n";
        // ret += ", " + format("w_unit_barrier", w_unit_barrier);
        // ret += ", " + format("w_global_rescale_param", w_global_rescale_param);
        // ret += ", " + format("w_fit", w_fit);

        ret += "\n\n" + format("CURR_STEP: ", total_step);
        ret += ", " + format("outer_step", outer_step);
        ret += ", " + format("w_outer_step", w_outer_step);

        ret += ", " + format("inner_iter_max_steps", inner_iter_max_steps);
        ret += ", " + format("grad_tol", grad_tol);
        ret += ", " + format("xTol", xTol);
        ret += ", " + format("fTol", fTol);
        ret += ", " + format("reg", reg) + "\n";

        return ret;
    }

    double getGlobalScale() { return w_global_rescale_param; }

   private:
    // the glbal scale which was applied to the current weights.
    double w_global_rescale_param = 1.;
};

}   // namespace MiNT3D