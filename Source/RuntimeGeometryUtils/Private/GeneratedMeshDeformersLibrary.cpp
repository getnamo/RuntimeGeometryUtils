
#include "GeneratedMeshDeformersLibrary.h"
#include "DynamicMesh/DynamicMesh3.h"
#include "FrameTypes.h"
#include "DynamicMesh/MeshNormals.h"
#include "Async/ParallelFor.h"
#include "GeneratedMesh.h"
#include "HeightmapDeformersLibrary.h"
#include "Interface_CollisionDataProviderCore.h"


using namespace UE::Geometry;

UGeneratedMesh* UGeneratedMeshDeformersLibrary::DeformMeshAxisSinWave1D(UGeneratedMesh* MeshObj, float Magnitude, float Frequency, float FrequencyShift, FVector AxisIn, FVector UpIn)
{
	FVector3d Axis(AxisIn), UpVector(UpIn);
	Axis.Normalize();
	UpVector.Normalize();

	if (MeshObj)
	{
		MeshObj->EditMeshInPlace([&](FDynamicMesh3& Mesh)
		{
			ParallelFor(Mesh.MaxVertexID(), [&](int32 vid)
			{
				if (Mesh.IsVertex(vid))
				{
					FVector3d Pos = Mesh.GetVertex(vid);
					double Dot = Pos.Dot(Axis);
					Dot += FrequencyShift;
					FVector3d NewPos = Pos + Magnitude * FMathd::Sin(Frequency * Dot) * UpVector;
					Mesh.SetVertex(vid, NewPos);
				}
			});
		});
	}

	return MeshObj;
}



UGeneratedMesh* UGeneratedMeshDeformersLibrary::DeformMeshMove(UGeneratedMesh* MeshObj, FVector FromWorldLocation /*= FVector(0, 0, 0)*/, FVector ToWorldLocation /*= FVector(0, 0, 0)*/, float Radius /*= 1*/, float Hardness /*= 1*/, float Magnitude /*= 1*/)
{
	if (MeshObj)
	{
		MeshObj->EditMeshInPlace([&](FDynamicMesh3& Mesh)
		{
			ParallelFor(Mesh.MaxVertexID(), [&](int32 vid)
			{
				if (Mesh.IsVertex(vid))
				{
					FVector3d Pos = Mesh.GetVertex(vid);
					float Distance = (Pos - FromWorldLocation).Length();
					
					if (Distance < Radius)
					{
						float Ratio = FMath::Clamp((Radius - Distance) / Radius, Hardness, 1.f);	//hardness makes a bottom clamp
						float Scale = Ratio * Magnitude;	//magnitude scale whole operation
						//float Displacement = (ToWorldLocation - FromWorldLocation).Size();
						FVector3d NewPos = Pos + ((ToWorldLocation - FromWorldLocation) * Scale);
						Mesh.SetVertex(vid, NewPos);
					}
					
				}
			});
		});
	}
	return MeshObj;
}

UGeneratedMesh* UGeneratedMeshDeformersLibrary::DeformMeshHeightmap(UGeneratedMesh* MeshObj, UTexture2D* HeightmapTexture, float ZScale/*= 1*/, FVector DirectionNormal/*= FVector(0, 0, 1)*/)
{
	//Lock texture for read and get some metadata
	int32 BytesPerPixel = 4;

	EPixelFormat Format = HeightmapTexture->GetPixelFormat();
	if (Format != PF_B8G8R8A8)
	{
		UE_LOG(LogTemp, Log, TEXT("NB: Pixel format is not PF_R8G8B8A8: %d"), Format);
		//Todo: detect format and adjust BytesPerPixel

		BytesPerPixel = 1;
	}

	const int32 TextureWidth = HeightmapTexture->GetPlatformData()->Mips[0].SizeX;
	const int32 TextureHeight= HeightmapTexture->GetPlatformData()->Mips[0].SizeY;
	const int32 DataLength = TextureHeight * TextureWidth * BytesPerPixel;

	void* TextureDataPointer = HeightmapTexture->GetPlatformData()->Mips[0].BulkData.Lock(LOCK_READ_ONLY);

	if (MeshObj)
	{
		MeshObj->EditMeshInPlace([&](FDynamicMesh3& Mesh)
		{
			//grab bounds for mapped sampling
			//TODO: generalize to projection vector
			FAxisAlignedBox3d Bounds = Mesh.GetBounds();
			float MeshWidth = Bounds.Width();
			float MeshHeight = Bounds.Height();
			FVector3d Center = Bounds.Center();

			//UE_LOG(LogTemp, Log, TEXT("MeshSize: %1.3f, %1.3f"), MeshWidth, MeshHeight);

			//TODO: optimize loop
			ParallelFor(Mesh.MaxVertexID(), [&](int32 vid)
			{
				if (Mesh.IsVertex(vid))
				{
					FVector3d Pos = Mesh.GetVertex(vid);

					//Normalize sample
					float XSample = (Pos.X + (MeshHeight / 2)) / MeshHeight * TextureHeight;
					float YSample = (Pos.Y + (MeshWidth / 2)) / MeshWidth * TextureWidth;

					//UE_LOG(LogTemp, Log, TEXT("Sample: %1.3f, %1.3f"), XSample, YSample);

					float UnscaledZ = UHeightmapDeformersLibrary::HeightAtPixel(XSample, YSample, TextureDataPointer, TextureHeight, TextureWidth, BytesPerPixel);

					//TODO: support vector deform

					//Apply new delta
					FVector3d NewPos = Pos;
					NewPos.Z = Pos.Z + UnscaledZ * ZScale;

					Mesh.SetVertex(vid, NewPos);
				}
			});
		});
	}

	//re-lock heightmap
	HeightmapTexture->GetPlatformData()->Mips[0].BulkData.Unlock();

	return MeshObj;
}

UGeneratedMesh* UGeneratedMeshDeformersLibrary::DeformMeshAxisSinWaveRadial(UGeneratedMesh* MeshObj, float Magnitude, float Frequency, float FrequencyShift, FVector AxisIn)
{
	FVector3d Axis(AxisIn);
	Axis.Normalize();
	FFrame3d AxisFrame(FVector3d::Zero(), Axis);

	if (MeshObj)
	{
		MeshObj->EditMeshInPlace([&](FDynamicMesh3& Mesh)
		{
			ParallelFor(Mesh.MaxVertexID(), [&](int32 vid)
			{
				if (Mesh.IsVertex(vid))
				{
					FVector3d Pos = Mesh.GetVertex(vid);
					double Dot = Pos.Dot(Axis);
					Dot += FrequencyShift;
					double Displacement = Magnitude * FMathd::Sin(Frequency * (Dot + FrequencyShift));
					FVector3d PlaneVec = AxisFrame.ToPlane(Pos);
					PlaneVec.Normalize();
					FVector3d NewPos = Pos + Displacement * PlaneVec;
					Mesh.SetVertex(vid, NewPos);
				}
			});
		});
	}

	return MeshObj;
}



UGeneratedMesh* UGeneratedMeshDeformersLibrary::DeformMeshPerlinNoiseNormal(UGeneratedMesh* MeshObj, float Magnitude, float Frequency, FVector FrequencyShift, int32 RandomSeed)
{
	if (MeshObj)
	{
		MeshObj->EditMeshInPlace([&](FDynamicMesh3& Mesh)
		{
			FMath::SRandInit(RandomSeed);
			const float RandomOffset = 10000.0f * FMath::SRand();
			FVector3d Offset(RandomOffset, RandomOffset, RandomOffset);
			Offset += (FVector3d)FrequencyShift;

			FMeshNormals Normals(&Mesh);
			Normals.ComputeVertexNormals();

			ParallelFor(Mesh.MaxVertexID(), [&](int32 vid)
			{
				if (Mesh.IsVertex(vid))
				{
					FVector3d Pos = Mesh.GetVertex(vid);
					FVector NoisePos = (FVector)( (double)Frequency * (Pos + Offset) );
					float Displacement = Magnitude * FMath::PerlinNoise3D(Frequency * NoisePos);
					FVector3d NewPos = Pos + Displacement * Normals[vid];
					Mesh.SetVertex(vid, NewPos);
				}
			});
		});
	}

	return MeshObj;
}

UGeneratedMesh* UGeneratedMeshDeformersLibrary::DeformMeshPerlinAlongAxisNormal(UGeneratedMesh* MeshObj, float Magnitude /*= 1*/, float Frequency /*= 1*/, FVector FrequencyShift /*= FVector(0, 0, 0)*/, int32 RandomSeed /*= 31337*/, FVector DirectionNormal /*= FVector(0, 0, 1)*/)
{
	if (MeshObj)
	{
		MeshObj->EditMeshInPlace([&](FDynamicMesh3& Mesh)
			{
				FMath::SRandInit(RandomSeed);
				const float RandomOffset = 10000.0f * FMath::SRand();
				FVector3d Offset(RandomOffset, RandomOffset, RandomOffset);
				Offset += (FVector3d)FrequencyShift;

				FMeshNormals Normals(&Mesh);
				Normals.ComputeVertexNormals();

				ParallelFor(Mesh.MaxVertexID(), [&](int32 vid)
					{
						if (Mesh.IsVertex(vid))
						{
							FVector3d Pos = Mesh.GetVertex(vid);
							FVector NoisePos = (FVector)((double)Frequency * (Pos + Offset));
							float Displacement = Magnitude * FMath::PerlinNoise3D(Frequency * NoisePos);
							FVector VertexNormal = (FVector)Normals[vid];

							if (acosf(VertexNormal | DirectionNormal) < (PI/4))
							{
								FVector3d NewPos = Pos + Displacement * VertexNormal;
								Mesh.SetVertex(vid, NewPos);
							}							
						}
					});
			});
	}

	return MeshObj;
}


UGeneratedMesh* UGeneratedMeshDeformersLibrary::DeformMeshPerlinNormalPastHeight(UGeneratedMesh* MeshObj, float Magnitude /*= 1*/, float Frequency /*= 1*/, FVector FrequencyShift /*= FVector(0, 0, 0)*/, int32 RandomSeed /*= 31337*/, FVector HeightVector /*= FVector(0, 0, 1)*/)
{
	if (MeshObj)
	{
		MeshObj->EditMeshInPlace([&](FDynamicMesh3& Mesh)
			{
				FMath::SRandInit(RandomSeed);
				const float RandomOffset = 10000.0f * FMath::SRand();
				FVector3d Offset(RandomOffset, RandomOffset, RandomOffset);
				Offset += (FVector3d)FrequencyShift;

				FMeshNormals Normals(&Mesh);
				Normals.ComputeVertexNormals();

				ParallelFor(Mesh.MaxVertexID(), [&](int32 vid)
					{
						if (Mesh.IsVertex(vid))
						{
							FVector3d Pos = Mesh.GetVertex(vid);
							FVector NoisePos = (FVector)((double)Frequency * (Pos + Offset));
							float Displacement = Magnitude * FMath::PerlinNoise3D(Frequency * NoisePos);
							FVector VertexNormal = (FVector)Normals[vid];

							if (Pos.Z >= HeightVector.Z)
							{
								FVector3d NewPos = Pos + Displacement * FVector(0,0,1);
								Mesh.SetVertex(vid, NewPos);
							}
						}
					});
			});
	}

	return MeshObj;
}

UGeneratedMesh* UGeneratedMeshDeformersLibrary::SmoothMeshUniform(UGeneratedMesh* MeshObj, float Alpha, int32 Iterations)
{
	Alpha = FMathf::Clamp(Alpha, 0.0f, 1.0f);
	Iterations = FMath::Clamp(Iterations, 0, 100);

	if (MeshObj)
	{
		MeshObj->EditMeshInPlace([&](FDynamicMesh3& Mesh)
		{
			int32 NumV = Mesh.MaxVertexID();
			TArray<FVector3d> SmoothPositions;
			SmoothPositions.SetNum(NumV);

			for (int32 k = 0; k < Iterations; ++k)
			{
				ParallelFor(NumV, [&](int32 vid)
				{
					if (Mesh.IsVertex(vid))
					{
						FVector3d Centroid;
						Mesh.GetVtxOneRingCentroid(vid, Centroid);
						SmoothPositions[vid] = FMath::Lerp<FVector3d>(Mesh.GetVertex(vid), Centroid, (double)Alpha);
					}
				});
				for (int32 vid = 0; vid < NumV; ++vid)
				{
					if (Mesh.IsVertex(vid))
					{
						Mesh.SetVertex(vid, SmoothPositions[vid]);
					}
				}
			}
		});
	}

	return MeshObj;
}

