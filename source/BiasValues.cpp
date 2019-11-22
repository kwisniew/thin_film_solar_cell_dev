#include "../include/BiasValues.hpp"

template<int dim>
void
Built_In_Bias<dim>::
set_value(const double & bias_value)
{
	built_in_bias = bias_value;
}

template<int dim>
void
Built_In_Bias<dim>::
set_location(const double & bias_location)
{
	built_in_bias_location = bias_location;
}

template <int dim>
double 
Built_In_Bias<dim>::
value(const dealii::Point<dim> &p, 
	  const unsigned int ) const
{
	// potential applied at y = 0, a positive potential is a foward bias,
	// a negative potential is a reverse bias
	if(p[0] > (built_in_bias_location - 1e-10)/*p[0] == 0.0*/) //|| (p[1] == 0.0))
	{
//		std::cout << "Built in" << std::endl;
		return built_in_bias; // volts

	}
	else 
		return 0;
//	double return_value = 0.0;
//	for(unsigned int i = 0; i < dim; i++)
//		return_value += 20.0*p[i];
//	return return_value;
}

template<int dim>
void
Schottky_Bias<dim>::
set_location(const double & bias_location)
{
	Schottky_location = bias_location;
}

template<int dim>
void
Schottky_Bias<dim>::
set_value(const double & bias_value)
{
	Schottky_bias = bias_value;
}

template <int dim>
double 
Schottky_Bias<dim>::
value(const dealii::Point<dim> &p, 
	  const unsigned int ) const
{
	// potential applied at x = 0, a positive potential is a foward bias,
	// a negative potential is a reverse bias
	if(p[0] < (Schottky_location + 1e-10) )
	{
//		std::cout << "on schottky" << std::endl;
		return Schottky_bias; // volts

	}
	else 
		return 0;
}

template<int dim>
void
Applied_Bias<dim>::
set_value(const double & bias_value)
{
	applied_bias = bias_value;
}

template<int dim>
void
Applied_Bias<dim>::
set_location(const double & bias_location)
{
	applied_bias_location = bias_location;
}


template <int dim>
double 
Applied_Bias<dim>::
value(const dealii::Point<dim> &p, 
	  const unsigned int ) const
{
	// potential applied at y = 0, a positive potential is a foward bias,
	// a negative potential is a reverse bias
	if(p[0] > (applied_bias_location - 1e-10)/*p[0] == 0.0*/) //|| (p[1] == 0))
	{
//		std::cout << "Applied" << std::endl;
		return applied_bias; // volts
	}
	else 
		return 0;
//	double return_value = 0.0;
//	for(unsigned int i = 0; i < dim; i++)
//		return_value += 20.0*p[i];
//	return return_value;
}
