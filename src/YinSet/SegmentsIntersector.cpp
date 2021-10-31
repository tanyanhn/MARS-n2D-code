#include "SegmentsIntersector.h"

void SegmentsIntersector::eventEnqueue(const rVec &p, cpSeg cit, int type)
{
  auto ret = event.insert(make_pair(p, std::map<cpSeg,int>()));
  std::map<cpSeg,int> &v = ret.first->second;
  v.insert(make_pair(cit,type));
}

void SegmentsIntersector::findNewEvents(const Bundle *s1, const Bundle *s2)
{
  for(auto seg1 = s1->cbegin(); seg1!=s1->cend(); ++seg1)
    for(auto seg2 = s2->cbegin(); seg2!=s2->cend(); ++seg2)
      findNewEvents(*seg1,*seg2);
}

void SegmentsIntersector::findNewEvents(cpSeg seg1, cpSeg seg2)
{
  rVec i1,i2;
  int n = intersect(*seg1, *seg2,i1,i2,tol);
  if(n>=1 && pt_cmp.compare(eventPoint,i1) != 1) {
    eventEnqueue(i1,seg1,INTERIOR);
    eventEnqueue(i1,seg2,INTERIOR);
  }
  if(n==2 && pt_cmp.compare(eventPoint,i2) != 1) {
    eventEnqueue(i2,seg1,INTERIOR);
    eventEnqueue(i2,seg2,INTERIOR);
  }
}

void SegmentsIntersector::collectCp(cpSeg testor)
{
  StatusType::iterator i = addToStatus(testor);
  StatusType::iterator j;
  if(i != status.begin()) {
    j = i; j--;
    findNewEvents(*i,*j);
  }
  j = i; j++;
  if(j != status.end()) {
    findNewEvents(*i,*j);
  }
  delFromStatus(testor);
}

SegmentsIntersector::StatusType::iterator
SegmentsIntersector::addToStatus(cpSeg cit)
{
  Bundle *v = new Bundle;
  v->push_back(cit);
  std::pair<StatusType::iterator,bool> ret = status.insert(v);
  if(!ret.second) {  // an already existed bundle implies overlap
    delete v;
    // test against others in the bundle
    Bundle *bdl = *ret.first;
    for(auto i=bdl->cbegin(); i!=bdl->cend(); i++)
      findNewEvents(cit, *i);
    bdl->push_back(cit);
  }
  // save the reference to the status structure, will be useful in delFromStatus()
  backref[cit-lineSetBegin] = ret.first;
  return ret.first;
}

SegmentsIntersector::StatusType::iterator
SegmentsIntersector::delFromStatus(cpSeg cit)
{
  StatusType::iterator ret = backref[cit-lineSetBegin];
  StatusType::iterator r2r = ret; r2r++;
  Bundle *bdl = *ret;
  bdl->erase(std::find(bdl->begin(), bdl->end(), cit));
  if(bdl->empty()) {
    delete bdl;
    status.erase(ret);
  }
  return r2r;
}

auto SegmentsIntersector::findIntersections(const std::vector<Seg> &lineSet) -> ResultType
{
  ResultType result(pt_cmp);
  
  // initialization of the event queue
  for(cpSeg i = lineSet.begin(); i != lineSet.end(); i++) {
    int s = pt_cmp.compare(i->p[0],i->p[1]);
    eventEnqueue(i->p[0], i, s);
    eventEnqueue(i->p[1], i, -s);
  }

  while(!event.empty()) {
    // pop the top event
    QueueType::iterator topEvent = event.begin();
    eventPoint = topEvent->first; // inform the compare object
    const std::map<cpSeg,int> &incident = topEvent->second;

    // check if L(p) V C(p) is empty
    bool emptyLC = std::all_of(incident.begin(), incident.end(),
			       [](const std::pair<cpSeg,int> &x) -> bool
			       { return x.second == UPPER; });
    // if yes, recover the true C(p)
    if(emptyLC) {
      cpSeg i = incident.cbegin()->first; // an arbitrary segment from U(p)
      collectCp(i);
    }

    // remove L(p) V C(p)
    StatusType::iterator r2r = status.end();
    for(auto i=incident.cbegin(); i!=incident.cend(); i++) {
      if(i->second != UPPER) 	// except U(p)
        r2r = delFromStatus(i->first);
    }
    
    StatusType::iterator leftmost = status.end();
    StatusType::iterator rightmost = status.end();
    
    // insert C(p) V U(p)
    for(auto i=incident.cbegin(); i!=incident.cend(); i++) {
      if(i->second != LOWER) {	// except L(p)
        StatusType::iterator j = addToStatus(i->first);
        // keep track of the leftmost and the rightmost bundle among C(p) V U(p)
        if(leftmost == status.end() || bdl_cmp(*j,*leftmost))
          leftmost = j;
        if(rightmost == status.end() || bdl_cmp(*rightmost,*j))
          rightmost = j;
      }
    }

    // test bundles that become newly adjacent for intersections
    if(leftmost == status.end()) { // the case of U(p) V C(p) empty
      if(r2r != status.end() && r2r != status.begin()) {
        StatusType::iterator l2l = r2r; l2l--;
        findNewEvents(*l2l,*r2r);
      }
    } else {
      // compare the leftmost one among C(p) V U(p) to the one left to it
      StatusType::iterator j;
      if(leftmost != status.begin()) {
        j = leftmost; j--;
        findNewEvents(*j,*leftmost);
      }
      // compare the rightmost one ... to the one right to it
      j = rightmost; j++;
      if(j != status.end()) {
        findNewEvents(*j,*rightmost);
      }
    }

    // pour to the output
    if(incident.size() > 1) {
      auto i = result.insert(make_pair(eventPoint, std::vector<int>()));
      std::vector<int> &outInc = i.first->second;
      for(auto i=incident.cbegin(); i!=incident.cend(); i++) {
        int index = i->first - lineSet.cbegin();
        outInc.push_back(index);
      }
    }

    event.erase(topEvent);
  } // end while !event.empty()
  backref.clear();
  return result;
}
