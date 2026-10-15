/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/

// upstream CubeCollisionModel.cpp:46; 7c18e95d5c5f2839079892c69e7d89a313c79603
void CubeCollisionModel::resize(sofa::Size size)
{
    const auto size0 = this->size;
    if (size == size0) return;
    // reset parent
    CollisionModel* parent = getPrevious();
    while(parent != nullptr)
    {
        parent->resize(0);
        parent = parent->getPrevious();
    }
    this->core::CollisionModel::resize(size);
    this->elems.resize(size);
    this->parentOf.resize(size);
    // set additional indices
    for (sofa::Size i=size0; i<size; ++i)
    {
        this->elems[i].children.first=core::CollisionElementIterator(getNext(), i);
        this->elems[i].children.second=core::CollisionElementIterator(getNext(), i+1);
        this->parentOf[i] = i;
    }
}

// upstream CubeCollisionModel.cpp:69; 7c18e95d5c5f2839079892c69e7d89a313c79603
void CubeCollisionModel::setParentOf(sofa::Index childIndex, const Vec3& min, const Vec3& max)
{
    const sofa::Index i = parentOf[childIndex];
    elems[i].minBBox = min;
    elems[i].maxBBox = max;
    elems[i].coneAngle = 2*M_PI;
}

// upstream CubeCollisionModel.cpp:77; 7c18e95d5c5f2839079892c69e7d89a313c79603
void CubeCollisionModel::setParentOf(sofa::Index childIndex, const Vec3& min, const Vec3& max, const Vec3& normal, const SReal angle)
{
    const sofa::Index i = parentOf[childIndex];
    elems[i].minBBox = min;
    elems[i].maxBBox = max;

    elems[i].coneAxis = normal;
    elems[i].coneAngle = angle;
}

// upstream CubeCollisionModel.cpp:103; 7c18e95d5c5f2839079892c69e7d89a313c79603
Index CubeCollisionModel::addCube(Cube subcellsBegin, Cube subcellsEnd)
{
    const sofa::Index index = size;

    this->core::CollisionModel::resize(index + 1);
    elems.resize(index + 1);
    
    elems[index].subcells.first = subcellsBegin;
    elems[index].subcells.second = subcellsEnd;
    elems[index].children.first = core::CollisionElementIterator();
    elems[index].children.second = core::CollisionElementIterator();
    updateCube(index);
    return index;
}

// upstream CubeCollisionModel.cpp:118; 7c18e95d5c5f2839079892c69e7d89a313c79603
void CubeCollisionModel::updateCube(sofa::Index index)
{
    const std::pair<Cube,Cube>& subcells = elems[index].subcells;
    if (subcells.first != subcells.second)
    {
        Cube c = subcells.first;
        Vec3 minBBox = c.minVect();
        Vec3 maxBBox = c.maxVect();

        elems[index].coneAxis = c.getConeAxis();
        elems[index].coneAngle = c.getConeAngle();

        ++c;
        while(c != subcells.second)
        {
            const Vec3& cmin = c.minVect();
            const Vec3& cmax = c.maxVect();

            const SReal alpha = std::max<SReal>(elems[index].coneAngle, c.getConeAngle());
            if(alpha <= M_PI/2)
            {
                const SReal beta = acos(c.getConeAxis() *  elems[index].coneAxis);
                elems[index].coneAxis = (c.getConeAxis() + elems[index].coneAxis).normalized();
                elems[index].coneAngle = beta/2 + alpha;
            }
            else
                elems[index].coneAngle = 2*M_PI;
            
            for (int j=0; j<3; j++)
            {
                if (cmax[j] > maxBBox[j]) maxBBox[j] = cmax[j];
                if (cmin[j] < minBBox[j]) minBBox[j] = cmin[j];
            }
            ++c;
        }
        elems[index].minBBox = minBBox;
        elems[index].maxBBox = maxBBox;
    }
}

// upstream CubeCollisionModel.cpp:158; 7c18e95d5c5f2839079892c69e7d89a313c79603
void CubeCollisionModel::updateCubes()
{
    for (sofa::Index i=0; i<size; i++)
        updateCube(i);
}

// upstream CubeCollisionModel.cpp:214; 7c18e95d5c5f2839079892c69e7d89a313c79603
std::pair<core::CollisionElementIterator,core::CollisionElementIterator> CubeCollisionModel::getInternalChildren(sofa::Index index) const
{
    return elems[index].subcells;
}

// upstream CubeCollisionModel.cpp:219; 7c18e95d5c5f2839079892c69e7d89a313c79603
std::pair<core::CollisionElementIterator,core::CollisionElementIterator> CubeCollisionModel::getExternalChildren(sofa::Index index) const
{
    return elems[index].children;
}

// upstream CubeCollisionModel.cpp:229; 7c18e95d5c5f2839079892c69e7d89a313c79603
void CubeCollisionModel::computeBoundingTree(int maxDepth)
{

    dmsg_info() << ">CubeCollisionModel::computeBoundingTree(" << maxDepth << ")";
    std::list<CubeCollisionModel*> levels;
    levels.push_front(createPrevious<CubeCollisionModel>());
    for (int i=0; i<maxDepth; i++)
        levels.push_front(levels.front()->createPrevious<CubeCollisionModel>());
    CubeCollisionModel* root = levels.front();

    if (root->empty() || root->getPrevious() != nullptr)
    {
        // Tree must be reconstructed
        dmsg_info() << "Building Tree with depth " << maxDepth << " from " << size << " elements.";
        // First remove extra levels
        while(root->getPrevious()!=nullptr)
        {
            core::CollisionModel::SPtr m = root->getPrevious();
            root->setPrevious(m->getPrevious());
            if (m->getMaster()) m->getMaster()->removeSlave(m);
            //delete m;
            m.reset();
        }

        // Then clear all existing levels
        {
            for (const auto & level : levels)
                level->resize(0);
        }

        // Then build root cell
        dmsg_info() << "CubeCollisionModel: add root cube";
        root->addCube(Cube(this,0),Cube(this,size));
        // Construct tree by splitting cells along their biggest dimension
        auto it = levels.begin();
        CubeCollisionModel* level = *it;
        ++it;
        int lvl = 0;
        while(it != levels.end())
        {
            dmsg_info() << "CubeCollisionModel: split level " << lvl;
            CubeCollisionModel* clevel = *it;
            clevel->elems.reserve(level->size*2);
            for(Cube cell = Cube(level->begin()); level->end() != cell; ++cell)
            {
                const std::pair<Cube,Cube>& subcells = cell.subcells();
                const sofa::Index ncells = subcells.second.getIndex() - subcells.first.getIndex();
                dmsg_info() << "CubeCollisionModel: level " << lvl << " cell " << cell.getIndex() << ": current subcells " << subcells.first.getIndex() << " - " << subcells.second.getIndex();
                if (ncells > 4)
                {
                    // Only split cells with more than 4 childs
                    // Find the biggest dimension
                    int splitAxis;
                    Vec3 l = cell.maxVect()-cell.minVect();
                    const sofa::Index middle = subcells.first.getIndex()+(ncells+1)/2;
                    if(l[0]>l[1])
                        if (l[0]>l[2])
                            splitAxis = 0;
                        else
                            splitAxis = 2;
                    else if (l[1]>l[2])
                        splitAxis = 1;
                    else
                        splitAxis = 2;

                    // Separate cells on each side of the median cell
                    const CubeSortPredicate sortpred(splitAxis);
                    std::stable_sort(elems.begin() + subcells.first.getIndex(), elems.begin() + subcells.second.getIndex(), sortpred);

                    // Create the two new subcells
                    const Cube cmiddle(this, middle);
                    sofa::Index c1 = clevel->addCube(subcells.first, cmiddle);
                    sofa::Index c2 = clevel->addCube(cmiddle, subcells.second);
                    dmsg_info() << "L" << lvl << " cell " << cell.getIndex() << " split along " << (splitAxis == 0 ? 'X' : splitAxis == 1 ? 'Y' : 'Z') << " in cell " << c1 << " size " << middle - subcells.first.getIndex() << " and cell " << c2 << " size " << subcells.second.getIndex() - middle << ".";
                    //level->elems[cell.getIndex()].subcells = std::make_pair(Cube(clevel,c1),Cube(clevel,c2+1));
                    level->elems[cell.getIndex()].subcells.first = Cube(clevel,c1);
                    level->elems[cell.getIndex()].subcells.second = Cube(clevel,c2+1);
                }
            }
            ++it;
            level = clevel;
            ++lvl;
        }
        if (!parentOf.empty())
        {
            // Finally update parentOf to reflect new cell order
            for (sofa::Size i=0; i<size; i++)
                parentOf[elems[i].children.first.getIndex()] = i;
        }
    }
    else
    {
        // Simply update the existing tree, starting from the bottom
        int lvl = 0;
        for (auto it = levels.rbegin(); it != levels.rend(); ++it)
        {
            dmsg_info() << "CubeCollisionModel: update level " << lvl;
            (*it)->updateCubes();
            ++lvl;
        }
    }
    dmsg_info() << "<CubeCollisionModel::computeBoundingTree(" << maxDepth << ")";
}
/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/


// upstream BVHNarrowPhase.cpp:42; 7c18e95d5c5f2839079892c69e7d89a313c79603
bool BVHNarrowPhase::isSelfCollision(core::CollisionModel* cm1, core::CollisionModel* cm2)
{
    return (cm1->getContext() == cm2->getContext());
}

// upstream BVHNarrowPhase.cpp:47; 7c18e95d5c5f2839079892c69e7d89a313c79603
void BVHNarrowPhase::addCollisionPair(const std::pair<core::CollisionModel*, core::CollisionModel*>& cmPair)
{
    core::CollisionModel *cm1 = cmPair.first; //->getNext();
    core::CollisionModel *cm2 = cmPair.second; //->getNext();

    if (!cm1->isSimulated() && !cm2->isSimulated())
        return;

    if (cm1->empty() || cm2->empty())
        return;

    core::CollisionModel *finestCollisionModel1 = cm1->getLast();//get the finest CollisionModel which is not a CubeModel
    core::CollisionModel *finestCollisionModel2 = cm2->getLast();

    const bool selfCollision = isSelfCollision(finestCollisionModel1, finestCollisionModel2);

    bool swapModels = false;
    core::collision::ElementIntersector* finestIntersector = intersectionMethod->findIntersector(finestCollisionModel1, finestCollisionModel2, swapModels);//find the method for the finest CollisionModels
    if (finestIntersector == nullptr)
        return;
    if (swapModels)
    {
        std::swap(cm1, cm2);
        std::swap(finestCollisionModel1, finestCollisionModel2);
    }

    // NOTE: outputs is a reference to a pointer! The original pointer resides in a map in NarrowPhaseDetection and will be modified in beginIntersect
    sofa::core::collision::DetectionOutputVector*& outputs = this->getDetectionOutputs(finestCollisionModel1, finestCollisionModel2);

    finestIntersector->beginIntersect(finestCollisionModel1, finestCollisionModel2, outputs);//creates outputs if null

    if (finestCollisionModel1 == cm1 || finestCollisionModel2 == cm2)
    {
        // The last model also contains the root element -> it does not only contains the final level of the tree
        finestCollisionModel1 = nullptr;
        finestCollisionModel2 = nullptr;
        finestIntersector = nullptr;
    }

    // Queue used for the iterative form of a tree traversal, avoiding the recursive form
    std::queue< TestPair > externalCells;
    initializeExternalCells(cm1, cm2, externalCells);

    core::collision::ElementIntersector* intersector = nullptr;
    MirrorIntersector mirror;
    cm1 = nullptr; // force later init of intersector
    cm2 = nullptr;

    while (!externalCells.empty())
    {
        TestPair root = externalCells.front();
        externalCells.pop();

        processExternalCell(root,
                            cm1, cm2,
                            intersector,
                            {finestCollisionModel1, finestCollisionModel2, finestIntersector, selfCollision},
                            &mirror, externalCells, outputs);
    }
}

// upstream BVHNarrowPhase.cpp:108; 7c18e95d5c5f2839079892c69e7d89a313c79603
void BVHNarrowPhase::initializeExternalCells(
        core::CollisionModel *cm1,
        core::CollisionModel *cm2,
        std::queue<TestPair>& externalCells)
{
    //See CollisionModel::getInternalChildren(Index), CollisionModel::getExternalChildren(Index) and definition of CollisionModel class
    const CollisionIteratorRange internalChildren1 = cm1->begin().getInternalChildren();
    const CollisionIteratorRange internalChildren2 = cm2->begin().getInternalChildren();
    const CollisionIteratorRange externalChildren1 = cm1->begin().getExternalChildren();
    const CollisionIteratorRange externalChildren2 = cm2->begin().getExternalChildren();

    const auto addToExternalCells = [&externalCells](
            const CollisionIteratorRange& children1,
            const CollisionIteratorRange& children2)
    {
        if (!isRangeEmpty(children1) && !isRangeEmpty(children2))
        {
            externalCells.emplace(children1,children2);
        }
    };

    addToExternalCells(internalChildren1, internalChildren2);
    addToExternalCells(internalChildren1, externalChildren2);
    addToExternalCells(externalChildren1, internalChildren2);
    addToExternalCells(externalChildren1, externalChildren2);
}

// upstream BVHNarrowPhase.cpp:135; 7c18e95d5c5f2839079892c69e7d89a313c79603
void BVHNarrowPhase::processExternalCell(const TestPair &externalCell,
                                              core::CollisionModel *&cm1,
                                              core::CollisionModel *&cm2,
                                              core::collision::ElementIntersector *coarseIntersector,
                                              const FinestCollision &finest,
                                              MirrorIntersector *mirror,
                                              std::queue<TestPair> &externalCells,
                                              sofa::core::collision::DetectionOutputVector *&outputs) const
{
    const auto [collisionModel1, collisionModel2] = getCollisionModelsFromTestPair(externalCell);

    if (cm1 != collisionModel1 || cm2 != collisionModel2)//if the CollisionElements do not belong to cm1 and cm2, update cm1 and cm2
    {
        cm1 = collisionModel1;
        cm2 = collisionModel2;
        if (!cm1 || !cm2) return;

        bool swapModels = false;
        coarseIntersector = intersectionMethod->findIntersector(cm1, cm2, swapModels);

        if (coarseIntersector == nullptr)
        {
            msg_error() << "Unable to find coarseIntersector " << intersectionMethod->getName() << " for "<<cm1->getClassName()<<" - "<<cm2->getClassName();
        }

        if (swapModels)
        {
            mirror->intersector = coarseIntersector;
            coarseIntersector = mirror;
        }
    }

    if (coarseIntersector == nullptr)
        return;

    // Stack used for the iterative form of a tree traversal, avoiding the recursive form
    std::stack< TestPair > internalCells;
    internalCells.push(externalCell);

    while (!internalCells.empty())
    {
        TestPair current = internalCells.top();
        internalCells.pop();

        processInternalCell(current, coarseIntersector, finest, externalCells, internalCells, outputs, intersectionMethod);
    }
}

// upstream BVHNarrowPhase.cpp:183; 7c18e95d5c5f2839079892c69e7d89a313c79603
void BVHNarrowPhase::processInternalCell(const TestPair &internalCell,
                                         core::collision::ElementIntersector *coarseIntersector,
                                         const FinestCollision &finest,
                                         std::queue<TestPair> &externalCells,
                                         std::stack<TestPair> &internalCells,
                                         sofa::core::collision::DetectionOutputVector *&outputs,
                                         const sofa::core::collision::Intersection* currentIntersection)
{
    const auto [collisionModel1, collisionModel2] = getCollisionModelsFromTestPair(internalCell);

    if (collisionModel1 == finest.cm1 && collisionModel2 == finest.cm2) //the collision models are the finest ones
    {
        // Final collision pairs
        finalCollisionPairs(internalCell, finest.selfCollision, coarseIntersector, outputs, currentIntersection);
    }
    else
    {
        visitCollisionElements(internalCell, coarseIntersector, finest, externalCells, internalCells, outputs, currentIntersection);
    }
}

// upstream BVHNarrowPhase.cpp:204; 7c18e95d5c5f2839079892c69e7d89a313c79603
void BVHNarrowPhase::visitCollisionElements(const TestPair &root,
                                            core::collision::ElementIntersector *coarseIntersector,
                                            const FinestCollision &finest,
                                            std::queue<TestPair> &externalCells,
                                            std::stack<TestPair> &internalCells,
                                            sofa::core::collision::DetectionOutputVector *&outputs,
                                            const sofa::core::collision::Intersection* currentIntersection)
{
    const core::CollisionElementIterator begin1 = root.first.first;
    const core::CollisionElementIterator end1 = root.first.second;
    const core::CollisionElementIterator begin2 = root.second.first;
    const core::CollisionElementIterator end2 = root.second.second;
    
    for (auto it1 = begin1; it1 != end1; ++it1)
    {
        for (auto it2 = begin2; it2 != end2; ++it2)
        {
            if (coarseIntersector->canIntersect(it1, it2, currentIntersection))
            {
                // Need to test recursively
                // Note that an element cannot have both internal and external children

                TestPair newInternalTests(it1.getInternalChildren(), it2.getInternalChildren());

                if (!isRangeEmpty(newInternalTests.first))
                {
                    if (!isRangeEmpty(newInternalTests.second))
                    {
                        //both collision elements have internal children. They are added to the list
                        internalCells.push(std::move(newInternalTests));
                    }
                    else
                    {
                        //only the first collision element has internal children. The second collision element
                        //is kept as it is
                        newInternalTests.second = {it2, it2 + 1};
                        internalCells.push(std::move(newInternalTests));
                    }
                }
                else
                {
                    if (!isRangeEmpty(newInternalTests.second))
                    {
                        //only the second collision element has internal children. The first collision element
                        //is kept as it is
                        newInternalTests.first = {it1, it1 + 1};
                        internalCells.push(std::move(newInternalTests));
                    }
                    else
                    {
                        // end of both internal tree of elements.
                        // need to test external children
                        visitExternalChildren(it1, it2, coarseIntersector, finest, externalCells, outputs, currentIntersection);
                    }
                }
            }
        }
    }
}

// upstream BVHNarrowPhase.cpp:264; 7c18e95d5c5f2839079892c69e7d89a313c79603
void BVHNarrowPhase::visitExternalChildren(const core::CollisionElementIterator &it1,
                                           const core::CollisionElementIterator &it2,
                                           core::collision::ElementIntersector *coarseIntersector,
                                           const FinestCollision &finest,
                                           std::queue<TestPair> &externalCells,
                                           sofa::core::collision::DetectionOutputVector *&outputs,
                                           const sofa::core::collision::Intersection* currentIntersection)
{
    const TestPair externalChildren(it1.getExternalChildren(), it2.getExternalChildren());

    const bool isExtChildrenRangeEmpty1 = isRangeEmpty(externalChildren.first);
    const bool isExtChildrenRangeEmpty2 = isRangeEmpty(externalChildren.second);

    if (!isExtChildrenRangeEmpty1)
    {
        if (!isExtChildrenRangeEmpty2)
        {
            const auto [collisionModel1, collisionModel2] = getCollisionModelsFromTestPair(externalChildren);
            if (collisionModel1 == finest.cm1 && collisionModel2 == finest.cm2) //the collision models are the finest ones
            {
                finalCollisionPairs(externalChildren, finest.selfCollision, finest.intersector, outputs, currentIntersection);
            }
            else
            {
                externalCells.push(externalChildren);
            }
        }
        else
        {
            // only first element has external children
            // test them against the second element
            externalCells.emplace(externalChildren.first, std::make_pair(it2, it2 + 1));
        }
    }
    else if (!isExtChildrenRangeEmpty2)
    {
        // only second element has external children
        // test them against the first element
        externalCells.emplace(std::make_pair(it1, it1 + 1), externalChildren.second);
    }
    else
    {
        // No child -> final collision pair
        if (!finest.selfCollision || it1.canCollideWith(it2))
            coarseIntersector->intersect(it1, it2, outputs, currentIntersection);
    }
}

// upstream BVHNarrowPhase.cpp:312; 7c18e95d5c5f2839079892c69e7d89a313c79603
void BVHNarrowPhase::finalCollisionPairs(const TestPair& pair,
                                         bool selfCollision,
                                         core::collision::ElementIntersector* intersector,
                                         sofa::core::collision::DetectionOutputVector*& outputs,
                                         const sofa::core::collision::Intersection* currentIntersection)
{
    const core::CollisionElementIterator begin1 = pair.first.first;
    const core::CollisionElementIterator end1 = pair.first.second;
    const core::CollisionElementIterator begin2 = pair.second.first;
    const core::CollisionElementIterator end2 = pair.second.second;

    for (auto it1 = begin1; it1 != end1; ++it1)
    {
        for (auto it2 = begin2; it2 != end2; ++it2)
        {
            // Final collision pair
            if (!selfCollision || it1.canCollideWith(it2))
                intersector->intersect(it1, it2, outputs, currentIntersection);
        }
    }
}

// upstream BVHNarrowPhase.cpp:334; 7c18e95d5c5f2839079892c69e7d89a313c79603
std::pair<core::CollisionModel*, core::CollisionModel*> BVHNarrowPhase::getCollisionModelsFromTestPair(const TestPair& pair)
{
    auto* collisionModel1 = pair.first.first.getCollisionModel(); //get the first collision model
    auto* collisionModel2 = pair.second.first.getCollisionModel(); //get the second collision model
    return {collisionModel1, collisionModel2};
}

// upstream BVHNarrowPhase.cpp:341; 7c18e95d5c5f2839079892c69e7d89a313c79603
bool BVHNarrowPhase::isRangeEmpty(const CollisionIteratorRange& range)
{
    return range.first == range.second;
}
/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
bool LMDIntersector::testCubes(Cube& cube1, Cube& cube2, const core::collision::Intersection* currentIntersection)
{
    const auto& minVect1 = cube1.minVect();
    const auto& minVect2 = cube2.minVect();
    const auto& maxVect1 = cube1.maxVect();
    const auto& maxVect2 = cube2.maxVect();

    const auto alarmDist = currentIntersection->getAlarmDistance() + cube1.getContactDistance() + cube2.getContactDistance();

    for (int i = 0; i < 3; i++)
    {
        if (minVect1[i] > maxVect2[i] + alarmDist || minVect2[i] > maxVect1[i] + alarmDist)
            return false;
    }

    return true;
}
/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/

// upstream PointCollisionModel.inl:115
template<class DataTypes>
void PointCollisionModel<DataTypes>::computeBoundingTree(int maxDepth)
{
    CubeCollisionModel* cubeModel = createPrevious<CubeCollisionModel>();
    const auto npoints = mstate->getSize();
    bool updated = false;
    if (npoints != size)
    {
        resize(npoints);
        updated = true;
    }
    if (updated) cubeModel->resize(0);
    if (!isMoving() && !cubeModel->empty() && !updated) return; // No need to recompute BBox if immobile

    if (d_computeNormals.getValue()) updateNormals();

    cubeModel->resize(size);
    if (!empty())
    {
        //VecCoord& x =mstate->read(core::vec_id::read_access::position)->getValue();
        const SReal distance = this->d_contactDistance.getValue();
        for (sofa::Size i=0; i<size; i++)
        {
            TPoint<DataTypes> p(this,i);
            const type::Vec3& pt = p.p();
            cubeModel->setParentOf(i, pt - type::Vec3(distance,distance,distance), pt + type::Vec3(distance,distance,distance));
        }
        cubeModel->computeBoundingTree(maxDepth);
    }
}
/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/

// upstream LineCollisionModel.inl:461
template<class DataTypes>
void LineCollisionModel<DataTypes>::computeBoundingTree(int maxDepth)
{
    CubeCollisionModel* cubeModel = createPrevious<CubeCollisionModel>();
    updateFromTopology();
    if (needsUpdate) cubeModel->resize(0);
    if (!isMoving() && !cubeModel->empty() && !needsUpdate) return; // No need to recompute BBox if immobile

    needsUpdate = false;

    cubeModel->resize(size);
    if (!empty())
    {
        const SReal distance = this->d_contactDistance.getValue();
        const auto& positions = this->mstate->read(core::vec_id::read_access::position)->getValue();
        for (sofa::Size i=0; i<size; i++)
        {
            type::Vec3 minElem, maxElem;
            const type::Vec3& pt1 = positions[this->elems[i].p[0]];
            const type::Vec3& pt2 = positions[this->elems[i].p[1]];

            for (int c = 0; c < 3; c++)
            {
                minElem[c] = pt1[c];
                maxElem[c] = pt1[c];
                if (pt2[c] > maxElem[c]) maxElem[c] = pt2[c];
                else if (pt2[c] < minElem[c]) minElem[c] = pt2[c];
                minElem[c] -= distance;
                maxElem[c] += distance;
            }

            cubeModel->setParentOf(i, minElem, maxElem);
        }
        cubeModel->computeBoundingTree(maxDepth);
    }
}
/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/

// upstream TriangleCollisionModel.inl:217
template<class DataTypes>
void TriangleCollisionModel<DataTypes>::computeBoundingTree(int maxDepth)
{
    CubeCollisionModel* cubeModel = createPrevious<CubeCollisionModel>();

    // check first that topology didn't changed
    if (m_topology->getRevision() != m_topologyRevision)
        updateFromTopology();

    if (m_needsUpdate && !cubeModel->empty())
        cubeModel->resize(0);

    if (!isMoving() && !cubeModel->empty() && !m_needsUpdate)
        return; // No need to recompute BBox if immobile nor if mesh didn't change.

    // set to false to avoid excessive loop
    m_needsUpdate=false;

    type::Vec3 minElem, maxElem;
    const VecCoord& x = this->m_mstate->read(core::vec_id::read_access::position)->getValue();

    const bool calcNormals = d_computeNormals.getValue();

    cubeModel->resize(size);  // size = number of triangles
    if (!empty())
    {
        const SReal distance = (SReal)this->d_contactDistance.getValue();
        for (sofa::Size i=0; i<size; i++)
        {
            Element t(this,i);

            const type::Vec3& pt1 = x[t.p1Index()];
            const type::Vec3& pt2 = x[t.p2Index()];
            const type::Vec3& pt3 = x[t.p3Index()];

            for (int c = 0; c < 3; c++)
            {
                minElem[c] = pt1[c];
                maxElem[c] = pt1[c];
                if (pt2[c] > maxElem[c]) maxElem[c] = pt2[c];
                else if (pt2[c] < minElem[c]) minElem[c] = pt2[c];
                if (pt3[c] > maxElem[c]) maxElem[c] = pt3[c];
                else if (pt3[c] < minElem[c]) minElem[c] = pt3[c];
                minElem[c] -= distance;
                maxElem[c] += distance;
            }
            if (calcNormals)
            {
                // Also recompute normal vector
                t.n() = cross(pt2-pt1,pt3-pt1);
                t.n().normalize();
            }

            if(d_useCurvature.getValue())
                cubeModel->setParentOf(i, minElem, maxElem, t.n()); // define the bounding box of the current triangle
            else
                cubeModel->setParentOf(i, minElem, maxElem);
        }
        cubeModel->computeBoundingTree(maxDepth);
    }
}
