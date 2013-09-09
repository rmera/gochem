/*
 * puckering.go, part of gochem.
 *
 *
 * Copyright 2012 Raul Mera <rmera{at}chemDOThelsinkiDOTfi>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2.1 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General
 * Public License along with this program.  If not, see
 * <http://www.gnu.org/licenses/>.
 *
 *
 * Gochem was started at the laboratory for instruction in Swedish, Department of Chemistry,
 * University of Helsinki, Finland. It is currently developed at AK Ochsenfeld, LMU-Munich,
 * Germany.
 *
 *
 */
/***Dedicated to the long life of the Ven. Khenpo Phuntzok Tenzin Rinpoche***/

//Generate (and in the future, measure) different puckering conformations for 
//5-member rings.

package chem

import (
	"fmt"
	"math"
)

//generate a new conformation for a 5-ring from the carbons ring and
//the phi and q parameters. 
//Based on Cremer and Pople, J Am Chem Soc, 96, 1354, (1975).
func  Move5ring(phi, q float64, carbons *CoordMatrix) *CoordMatrix{
	const (
		N = 5.0
	)
	centeredold:=Zeros(N,3)
	centered,displacement,_:=MassCentrate(carbons,carbons,nil)
	fmt.Println(centered)
	Rp:=Zeros(1,3)
	Rpp:=Zeros(1,3)
	tmp:=Zeros(1,3)
	tmp2:=Zeros(1,3)
	for j:=0;j<3;j++{
		tmp.Set(0,0,centered.At(j,0))
		tmp.Set(0,1,centered.At(j,1))
		tmp.Set(0,2,centered.At(j,2))
		tmp2.Clone(tmp)
		F:=math.Sin(2*math.Pi*(float64(j))/N)
		F2:=math.Cos(2*math.Pi*(float64(j))/N)
		tmp.Scale(F,tmp)
		Rp.Add(Rp,tmp)
		tmp2.Scale(F2,tmp2)
		Rpp.Add(Rpp,tmp2)
	}
	fmt.Println(centered)
	Normal:=Cross3DRow(Rp,Rpp)
	Normal.Unit(Normal)
	ZTrans:=GetSwitchZ(Normal)
	centered.Mul(centered,ZTrans)
	oldZ:=NewCoords([]float64{0,0,1},1,3)
	oldY:=NewCoords([]float64{0,1,0},1,3)
	oldZ.Mul(oldZ,ZTrans)
	oldY.Mul(oldY,ZTrans)
	fmt.Println(centered,Normal,Rp,Rpp)
	first:=EmptyCoords()
	first.RowView(centered,0)
	YTrans,_:=GetRotateToNewY(first)
	centered.Mul(centered,YTrans)
	oldY.Mul(oldY,YTrans)
	fmt.Println(centered)
	centeredold.Clone(centered)
	//Now we have the molecule correctly oriented.
	for j:=0;j<3;j++{
		z:=math.Sqrt(2/N)*q*math.Cos(phi + 4*math.Pi*float64(j)/N)
		centered.Set(j,2,z)
	}
	scaleXY(centeredold,centered,N)
	ZTrans=GetSwitchZ(oldZ)
	centered.Mul(centered,ZTrans)
	oldY.Mul(oldY,ZTrans)    //I think this and the previous  oldY.Mul(oldY,ZTrans) transformtions could be avoided
	YTrans,_=GetRotateToNewY(oldY)
	centered.Mul(centered,YTrans)
	centered.AddRow(centered,displacement)
	return centered
}

func scaleXY(old, curr *CoordMatrix, N int){
	oldr:=EmptyCoords()
	newr:=EmptyCoords()
	fmt.Println("old", curr)
	tmp:=Zeros(1,3)
	for i:=0;i<N;i++{
		oldr.RowView(old,i)
		newr.RowView(curr,i)
		size:=math.Sqrt(math.Pow(oldr.Norm(2),2)-math.Pow(newr.At(0,2),2))
		tmp.Clone(newr)
		tmp.Set(0,2,0.0)
		fmt.Println("OldNorms",oldr.Norm(2),newr.Norm(2))
		tmp.Scale(size/tmp.Norm(2),tmp)
		newr.Set(0,0,tmp.At(0,0))
		newr.Set(0,1,tmp.At(0,1))
		fmt.Println("Norms!", oldr.Norm(2),newr.Norm(2))
//		diff:=(math.Abs(newr.At(0,2))-math.Abs(oldr.At(0,2)))/2   //(newr.Norm(2)-oldr.Norm(2))/2.0
//		fmt.Println("diff", diff)
//		tmp.Set(0,0,diff)
//		tmp.Set(0,1,diff)
//		newr.Sub(newr,tmp)
	}
	fmt.Println("new",curr)

}

