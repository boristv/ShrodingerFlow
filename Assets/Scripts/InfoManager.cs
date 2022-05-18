using UnityEngine;
using UnityEngine.UI;

public class InfoManager : MonoBehaviour
{
    [SerializeField] private Text _iterationNumber;
    [SerializeField] private SFBase _sfBase;

    public void Update()
    {
        _iterationNumber.text = _sfBase.Iterator.ToString();
    }
}
