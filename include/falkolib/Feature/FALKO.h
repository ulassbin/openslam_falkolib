/**
 * FALKOLib - Fast Adaptive Laser Keypoint Orientation-invariant
 * Copyright (C) 2016 Fabjan Kallasi and Dario Lodi Rizzini.
 *
 * This file is part of FALKOLib.
 *
 * FALKOLib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * FALKOLib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with FALKOLib.  If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once

#include <falkolib/Feature/Keypoint.h>

namespace falkolib{
	
	/**
	 * @brief class representing a FALKO keypoint
	 * 
	 * FALKO keypoint extends simple keypoint interface.
	 * The index in the original scan, corner orientation and neighborhood search radius are also preserved.
	 */
	class FALKO : public Keypoint{
	public:
		int index;
		double radius;
		double orientation;
		double score;

		// Triangle variables
		double left_r; 
		int left_index;
		double right_r;
		int right_index;

		double left_x, left_y;
		double right_x, right_y;

		// Area variables
		double area;
		double base_len;
		double base_height;
	};
}