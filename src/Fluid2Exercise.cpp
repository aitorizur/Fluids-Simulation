#include "Scene.h"
#include "Numeric/PCGSolver.h"

namespace asa
{

namespace
{
////////////////////////////////////////////////
// Add any reusable classes or functions HERE //
////////////////////////////////////////////////


}  // namespace

void Fluid2::fluidAdvection(const float dt)
{
    // Ink Advection
    Array2<Vector3> inkRGBCopy(inkRGB);
    Index2 inkRGBSize = inkRGB.getSize();

    for (int i = 0; i < inkRGBSize.x; i++) 
    {
        for (int j = 0; j < inkRGBSize.y; j++) 
        {
            Index2 id = Index2(i, j);
            Vector2 cellPosition = grid.getCellPos(id);

            float iNext = floor(clamp(i + 1, 0, inkRGBSize.x - 1));
            float jNext = floor(clamp(j + 1, 0, inkRGBSize.y - 1));
            float vx = velocityX[id] + velocityX[Index2(iNext, j)] * 0.5f;
            float vy = velocityY[id] + velocityY[Index2(i, jNext)] * 0.5f;

            Vector2 prevPosition = cellPosition - dt * Vector2(vx, vy);
            Vector2 prevId = grid.getCellIndex(prevPosition);

            float clampPrevX = floor(clamp(prevId.x, 0, inkRGBSize.x - 1));
            float clampX = floor(clamp(prevId.x + 1, 0, inkRGBSize.x - 1));
            float clampPrevY = floor(clamp(prevId.y, 0, inkRGBSize.y - 1));
            float clampY = floor(clamp(prevId.y + 1, 0, inkRGBSize.y - 1));

            Vector3 aa = inkRGBCopy[Index2(clampPrevX, clampPrevY)];
            Vector3 ba = inkRGBCopy[Index2(clampX, clampPrevY)];
            Vector3 ab = inkRGBCopy[Index2(clampPrevX, clampY)];
            Vector3 bb = inkRGBCopy[Index2(clampX, clampY)];

            float t = prevId.x - floor(prevId.x);
            float s = prevId.y - floor(prevId.y);
            Vector3 bilerpInk = bilerp(aa, ba, ab, bb, t, s);
            inkRGB[id] = bilerpInk;
        }
    }

    // Velocity advection
    Array2<float> velocityXCopy(velocityX);
    Array2<float> velocityYCopy(velocityY);
    Index2 sizeVelocityX = velocityX.getSize();
    Index2 sizeVelocityY = velocityY.getSize();

    // Velocity X
    for (int i = 0; i < sizeVelocityX.x; i++) 
    {
        for (int j = 0; j < sizeVelocityX.y; j++) 
        {
            Index2 id = Index2(i, j);

            Vector2 facePosition = grid.getFacePosX(id);

            int clampI = clamp(i, 0, sizeVelocityY.x - 1);
            int clampIPrev = clamp(i - 1, 0, sizeVelocityY.x - 1);
            int clampJ = clamp(j, 0, sizeVelocityY.y - 1);
            int clampJPrev = clamp(j - 1, 0, sizeVelocityY.y - 1);

            float u = velocityXCopy[id];

            float v0 = velocityYCopy[Index2(clampIPrev, clampJ)];
            float v1 = velocityYCopy[Index2(clampI, clampJ)];
            float v2 = velocityYCopy[Index2(clampIPrev, clampJPrev)];
            float v3 = velocityYCopy[Index2(clampI, clampJPrev)];

            float v = (v0 + v1 + v2 + v3) * 0.25f;

            Vector2 prevPosition = facePosition - dt * Vector2(u, v);
            Vector2 prevId = grid.getCellIndex(prevPosition);

            float clampVIX = floor(clamp(prevId.x, 0, sizeVelocityX.x - 1));
            float clampVIX1 = floor(clamp(prevId.x + 1, 0, sizeVelocityX.x - 1));
            float clampVIY = floor(clamp(prevId.y, 0, sizeVelocityX.y - 1));
            float clampVIY1 = floor(clamp(prevId.y + 1, 0, sizeVelocityX.y - 1));

            float aa = velocityXCopy[Index2(clampVIX, clampVIY)];
            float ba = velocityXCopy[Index2(clampVIX1, clampVIY)];
            float ab = velocityXCopy[Index2(clampVIX, clampVIY1)];
            float bb = velocityXCopy[Index2(clampVIX1, clampVIY1)];

            float t = prevId.x - floor(prevId.x);
            float s = prevId.y - floor(prevId.y);

            velocityX[id] = bilerp(aa, ba, ab, bb, t, s);
        }
    }

    // Velocity Y
    for (int i = 0; i < sizeVelocityY.x; i++) 
    {
        for (int j = 0; j < sizeVelocityY.x; j++) 
        {
            Index2 id = Index2(i, j);

            Vector2 facePosition = grid.getFacePosY(id);

            int clampI = clamp(i, 0, sizeVelocityX.x - 1);
            int clampINext = clamp(i + 1, 0, sizeVelocityX.x - 1);
            int clampIPrev = clamp(i - 1, 0, sizeVelocityX.x - 1);
            int clampJ = clamp(j, 0, sizeVelocityX.y - 1);
            int clampJNext = clamp(j + 1, 0, sizeVelocityX.y - 1);
            int clampJPrev = clamp(j - 1, 0, sizeVelocityX.y - 1);

            float v = velocityYCopy[id];

            float u0 = velocityXCopy[Index2(clampI, clampJPrev)];
            float u1 = velocityXCopy[Index2(clampI, clampJ)];
            float u2 = velocityXCopy[Index2(clampINext, clampJPrev)];
            float u3 = velocityXCopy[Index2(clampINext, clampJ)];

            float u = (u0 + u1 + u2 + u3) * 0.25f;

            Vector2 prevPosition = facePosition - dt * Vector2(u, v);
            Vector2 prevId = grid.getCellIndex(prevPosition);

            float clampVIX = floor(clamp(prevId.x, 0, sizeVelocityY.x - 1));
            float clampVIX1 = floor(clamp(prevId.x + 1, 0, sizeVelocityY.x - 1));
            float clampVIY = floor(clamp(prevId.y, 0, sizeVelocityY.y - 1));
            float clampVIY1 = floor(clamp(prevId.y + 1, 0, sizeVelocityY.y - 1));

            float aa = velocityYCopy[Index2(clampVIX, clampVIY)];
            float ba = velocityYCopy[Index2(clampVIX1, clampVIY)];
            float ab = velocityYCopy[Index2(clampVIX, clampVIY1)];
            float bb = velocityYCopy[Index2(clampVIX1, clampVIY1)];

            float t = prevId.x - floor(prevId.x);
            float s = prevId.y - floor(prevId.y);

            velocityY[id] = bilerp(aa, ba, ab, bb, t, s);
        }
    }
}

void Fluid2::fluidEmission()
{
    if (Scene::testcase >= Scene::SMOKE) 
    {
        int startRedX = 39;
        int startGreenX = 44;
        int startBlueX = 54;
        int startY = 0;
        int endY = 15;
        int endX = 59;

        Vector3 color = Vector3(1.0f, 0.0f, 0.0f);

        for (int i = startRedX; i < endX; i++) 
        {
            if (i >= startGreenX) 
            {
                color = Vector3(0.0f, 1.0f, 0.0f);
            }
            if (i >= startBlueX) 
            {
                color = Vector3(0.0f, 0.0f, 1.0f);
            }

            for (int j = startY; j < endY; j++) 
            {
                // Give starting color and velocity
                Index2 ind(i, j);
                inkRGB[ind] = color;
                velocityY[ind] = 10.0f;
            }
        }
    }
}

void Fluid2::fluidVolumeForces(const float dt)
{
    if (Scene::testcase >= Scene::SMOKE) 
    {
        Index2 sizeY = velocityY.getSize();

        for (int i = 0; i < sizeY.x; ++i) 
        {
            for (int j = 0; j < sizeY.y; ++j) 
            {
                Index2 id(i, j);
                velocityY[id] += dt * Scene::kGravity;  // Add gravity
            }
        }
    }
}

void Fluid2::fluidViscosity(const float dt)
{
    if (Scene::testcase >= Scene::SMOKE) 
    {
        Array2<float> copyVelocityX(velocityX);
        Array2<float> copyVelocityY(velocityY);

        Index2 sizeX = velocityX.getSize();
        Index2 sizeY = velocityY.getSize();

        float squareDx = pow(grid.getDx().x, 2);
        float squareDy = pow(grid.getDx().y, 2);

        // viscosity X
        for (int i = 0; i < sizeX.x; ++i) 
        {
            for (int j = 0; j < sizeX.y; ++j) 
            {
                Index2 id(i, j);

                Index2 idINext(clamp(i + 1, 0, sizeX.x - 1), j);
                Index2 idIPrev(clamp(i - 1, 0, sizeX.x - 1), j);
                Index2 idJNext(i, clamp(j + 1, 0, sizeX.y - 1));
                Index2 idjPrev(i, clamp(j - 1, 0, sizeX.y - 1));

                velocityX[id] += dt * Scene::kViscosity / Scene::kDensity * ((copyVelocityX[idINext] - 2.0f * copyVelocityX[id] + copyVelocityX[idIPrev]) / squareDx +
                                        (copyVelocityX[idJNext] - 2.0f * copyVelocityX[id] + copyVelocityX[idjPrev]) / squareDy);
            }
        }

        // viscosity Y
        for (int i = 0; i < sizeY.x; ++i) 
        {
            for (int j = 0; j < sizeY.y; ++j) 
            {
                Index2 id(i, j);

                Index2 idINext(clamp(i - 1, 0, sizeY.x - 1), j);
                Index2 idIPrev(clamp(i + 1, 0, sizeY.x - 1), j);
                Index2 idJNext(i, clamp(j - 1, 0, sizeY.y - 1));
                Index2 idJPrev(i, clamp(j + 1, 0, sizeY.y - 1));

                velocityY[id] += dt * Scene::kViscosity / Scene::kDensity * ((copyVelocityY[idINext] - 2.0f * copyVelocityY[id] + copyVelocityY[idIPrev]) * squareDx +
                                        (copyVelocityY[idJNext] - 2.0f * copyVelocityY[id] + copyVelocityY[idJPrev]) * squareDy);
            }
        }
    }
}

void Fluid2::fluidPressureProjection(const float dt)
{
    if (Scene::testcase >= Scene::SMOKE) 
    {
        float dx = grid.getDx().x;
        float dy = grid.getDx().y;
        float inverseSquareDx = 1.0f / pow(grid.getDx().x, 2);
        float inverseSquareDy = 1.0f / pow(grid.getDx().y, 2);

        Index2 pressureSize = pressure.getSize();
        Index2 velocityXSize = velocityX.getSize();
        Index2 velocityYSize = velocityY.getSize();

        // Boundaries
        for (int j = 0; j < velocityXSize.y; ++j) 
        {
            velocityX[Index2(0, j)] = 0.0f;
            velocityX[Index2(velocityXSize.x - 1, j)] = 0.0f;
        }
        for (int i = 0; i < velocityYSize.x; ++i) 
        {
            velocityY[Index2(i, 0)] = 0.0f;
        }

        // Matrix b
        std::vector<double> b(pressureSize.x * pressureSize.y);
        for (int i = 0; i < pressureSize.x; ++i) 
        {
            for (int j = 0; j < pressureSize.y; ++j) 
            {
                Index2 id(i, j);
                b[pressure.getLinearIndex(i, j)] = -Scene::kDensity / dt * ((velocityX[Index2(i + 1, j)] - velocityX[id]) / dx +
                                                                           (velocityY[Index2(i, j + 1)] - velocityY[id]) / dy);
            }
        }

        // Matrix A
        SparseMatrix<double> A(pressureSize.x * pressureSize.y, 10);
        for (int i = 0; i < pressureSize.x; ++i) 
        {
            for (int j = 0; j < pressureSize.y; ++j) 
            {
                int id = pressure.getLinearIndex(i, j);
                if (i > 0) 
                {
                    int idIPrev = pressure.getLinearIndex(i - 1, j);
                    A.add_to_element(id, id, inverseSquareDx);
                    A.add_to_element(id, idIPrev, -1.0 * inverseSquareDx);
                }
                if (i < pressureSize.x - 1) 
                {
                    int idINext = pressure.getLinearIndex(i + 1, j);
                    A.add_to_element(id, id, inverseSquareDx);
                    A.add_to_element(id, idINext, -1.0 * inverseSquareDx);
                }
                if (j > 0) 
                {
                    int idJPrev = pressure.getLinearIndex(i, j - 1);
                    A.add_to_element(id, id, inverseSquareDy);
                    A.add_to_element(id, idJPrev, -1.0 * inverseSquareDy);
                }

                A.add_to_element(id, id, inverseSquareDy);
                if (j < pressureSize.y - 1) 
                {
                    int idJNext = pressure.getLinearIndex(i, j + 1);
                    A.add_to_element(id, idJNext, -1.0 * inverseSquareDy);
                }
            }
        }

        // Pressure solve
        PCGSolver<double> solver;
        solver.set_solver_parameters(0.0001, 100000);

        double residualOut;
        int iterationsOut;
        std::vector<double> result(pressureSize.x * pressureSize.y);
        solver.solve(A, b, result, residualOut, iterationsOut);

        // Pressure
        int n = pressureSize.x * pressureSize.y;
        for (int i = 0; i < n; ++i) 
        {
            pressure[i] = (float)result[i];
        }

        // Pressure gradient
        float constant = -dt / Scene::kDensity;

        for (int i = 1; i < velocityXSize.x - 1; ++i) 
        {
            for (int j = 0; j < velocityXSize.y; ++j) 
            {
                Index2 id(i, j);
                float pressureGrad = (pressure[id] - pressure[Index2(i - 1, j)]) / dx;
                velocityX[id] += constant * pressureGrad;
            }
        }

        for (int i = 0; i < velocityYSize.x; ++i) 
        {
            for (int j = 1; j < velocityYSize.y - 1; ++j) 
            {
                Index2 id(i, j);
                float pressureGrad = (pressure[id] - pressure[Index2(i, j - 1)]) / dy;
                velocityY[id] += constant * pressureGrad;
            }
        }
    }
}
}  // namespace asa
